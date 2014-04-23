#include "StdAfx.h"
#include "CycloDetector.h"

Stroke::Stroke()
{
	strokePoints.clear();

	direction.x=0;
	direction.y=0;

	speed.x=0;
	speed.y=0;

	startTick=0;
	endTick=0;

	frequence=0;
	gain=1;

}
Stroke::~Stroke()
{
}

vector<TouchPoint> Stroke::getStrokePoints()
{
	return this->strokePoints;
}

void Stroke::addStrokePoint(TouchPoint point)
{
	this->strokePoints.push_back(point);
}

Point2f Stroke::getStrokeDirection()
{
	return this->direction;
}

void Stroke::setStrokeDirection(Point2f direction)
{
	this->direction=direction;
}

Point2f Stroke::getStrokeSpeed()
{
	return this->speed;
}

void Stroke::setStrokeSpeed(Point2f speed)
{
	this->speed=speed;
}

void Stroke::setStartTick(int64 startTick)
{
	this->startTick=startTick;
}

void Stroke::setEndTick(int64 endTick)
{
	this->endTick=endTick;

	float Te=(this->endTick-this->startTick)/getTickFrequency();
	this->frequence=1/(2*Te);
	//std::cout<<"Te"<<Te<<std::endl;
	TouchPoint startPoint=this->strokePoints[0];
	TouchPoint endPoint=this->strokePoints[this->strokePoints.size()-1];

	this->speed.x=(endPoint.x-startPoint.x)/Te;
	this->speed.y=(endPoint.y-startPoint.y)/Te;
}

void Stroke::setFrequence(float frequence)
{
	this->frequence=frequence;
}

float Stroke::getFrequence()
{
	return this->frequence;
}

void Stroke::setGain(vector<Stroke*> strokes)
{
	float k=(float)2/3;
	
	int strokeCount=0;
	float frequenceSum=0;

	for(int i=strokes.size()-1;i>=0;--i)
	{
		strokeCount++;
		if(strokeCount<=2)
		frequenceSum+=strokes[i]->getFrequence();
		else
			break;
	}

	float panGain=k*frequenceSum/strokeCount;
	//std::cout<<"sum "<<panGain<<" strokeCount "<<strokeCount<<std::endl;
	1>panGain? gain=1:gain=k/strokeCount*frequenceSum;
}

float Stroke::getGain()
{
	return this->gain;
}

CycloDetector::CycloDetector(void)
{
	dragStart=false;
	firstDragEnd=false;
	dragDirection.x=0;
	dragDirection.y=0;
}

CycloDetector::~CycloDetector(void)
{
}

void CycloDetector::init()
{
	touch.Init();
}

void CycloDetector::detectTouch()
{ 
	Mat touchImage;
	touchImage=Mat::zeros(480,640,CV_8UC3);

	lastTouchPoints=touchPoints;

	if(touch.checkTouch()==true)
	{
		for(int i=0;i<touch.getTouchPointList().size();++i)
		{
			Point touchPoint(touch.getTouchPointList()[i].x*640/touch.getResolutionX(),
				touch.getTouchPointList()[i].y*480/touch.getResolutionY());

			circle(touchImage,touchPoint,3,Scalar(0,255,0));
		}

		touchPoints=touch.getTouchPointList();
	}
	else
		touchPoints.clear();

	imshow("touchImage",touchImage);
}

float CycloDetector::getDistance(TouchPoint point1,TouchPoint point2)
{
	float disX=(float)(point1.x-point2.x);
	float disY=(float)(point1.y-point2.y);

	float dis=disX*disX+disY*disY;
	dis=powf(dis,0.5);

	return dis;
}

Point2f CycloDetector::getDirection(TouchPoint endPoint,TouchPoint startPoint)
{
	Point2f direction(endPoint.x-startPoint.x,endPoint.y-startPoint.y);
	float dis=getDistance(endPoint,startPoint);
	if(dis!=0)
	{
		direction.x/=getDistance(endPoint,startPoint);
		direction.y/=getDistance(endPoint,startPoint);
	}

	return direction;
}

bool CycloDetector::checkReverse(Point2f direction,Point2f lastDirection)
{
	float dotValue=direction.dot(lastDirection);
	
	bool reverse;
	dotValue>=0? reverse=false:reverse=true;
	return reverse;
}

void CycloDetector::cycloPan()
{
	Mat testImage=Mat::zeros(480,640,CV_8UC3);
	
	//Cancel the cyclo pan if not only one finger is put on the screen
	if(touchPoints.size()!=1)
		cancelPan();
	else
	{
		//Get the tick count for each frame
		int64 tickCount=getTickCount();
		
		if(dragStart==false)
		{
			//Start pan
			startPan(tickCount);
		}
		else
		{
			lastDragDirection=dragDirection;
			//Check the drag direction between this point and the last point
			dragDirection=getDirection(touchPoints[0],lastTouchPoints[0]);
			dragAmplitude=getDistance(touchPoints[0],lastTouchPoints[0]);


			/*std::cout<<"dragAmplitude "<<dragAmplitude<<std::endl;*/

			//When the first drag is still active
			if(firstDragEnd==false&&dragDirection!=Point2f(0,0))
			{
				//If the distance between the start point and the touch point in this frame is
				//larger than that between the start point and the last touch point, refresh the direction
				//of the first stroke
				if(getDistance(startDragPoint,touchPoints[0])>getDistance(startDragPoint,lastTouchPoints[0]))
				{
					stroke->setStrokeDirection(dragDirection);
					firstDragDirection=dragDirection;
				}
				/*stroke->setStrokeDirection(dragDirection);
				firstDragDirection=dragDirection;
				firstDrag=true;*/
			}

			//Check if the break time is larger than the threshold
			if(dragAmplitude==0)
			{
				
				breakTime+=(tickCount-lastTickCount)/getTickFrequency();
				
				//If the break time is larger then 0.2s, cancel the drag
				if(breakTime>0.2)
					cancelPan();
				else
					lastTickCount=tickCount;
			}
			//If there is no break in this frame
			else
			{
				breakTime=0;
				lastTickCount=tickCount;
			}

			//If the dragging direction is not reversed
			if(checkReverse(dragDirection,stroke->getStrokeDirection())==false)
			{
				//Add this point to the active stroke
				stroke->addStrokePoint(touchPoints[0]);
				//std::cout<<"gain "<<stroke->getGain()<<std::endl;
				Point2f translation=firstDragDirection*stroke->getGain()*dragAmplitude;
				translation.x/=10.0f;
				translation.y/=10.0f;
				//std::cout<<"amplitude"<<translation<<std::endl;
				testPoint+=translation;
			}
			//If the dragging direction is revesed
			else
			{
				//Set the end time of the stroke
				stroke->setEndTick(tickCount);
				
				//Add the active stroke to the vector of strokes
				strokes.push_back(stroke);
				//std::cout<<"speed "<<norm(stroke->getStrokeSpeed())<<std::endl;

				//If this is the first stroke, recalculate its drag direction
				if(strokes.size()==1)
				{
					firstDragEnd=true;
					firstDragDirection=getDirection(touchPoints[0],stroke->getStrokePoints()[0]);
					
				}
				//Check if the cylcoPan can start
				if(strokes.size()==1&&norm(stroke->getStrokeSpeed())<50)		
					cancelPan();
				else
				{
					//Refresh the active stroke
					stroke=new Stroke();
					stroke->setStartTick(tickCount);
					stroke->setStrokeDirection(dragDirection);
					stroke->addStrokePoint(touchPoints[0]);
					stroke->setGain(strokes);
					std::cout<<"reverse"<<std::endl;
				}
			}	
		}

		//for(int i=0;i<strokes.size();++i)
		//{
		//	float speed=norm(strokes[i]->getStrokeSpeed());
		//	std::cout<<"stroke "<<i<<" point count "<<strokes[i]->getStrokePoints().size()<<" speed "<<speed<<" gain "<<strokes[i]->getGain()<<std::endl;
		//}

		Point drawPoint(testPoint.x*640/touch.getResolutionX(),testPoint.y*480/touch.getResolutionY());
	    
		circle(testImage,Point(320,240),3,Scalar(0,255,255));
		circle(testImage,drawPoint,3,Scalar(0,255,0));
		imshow("testImage",testImage);
	}
}


void CycloDetector::cancelPan()
{
	dragStart=false;
	firstDragEnd=false;

	lastDragDirection.x=0;
	lastDragDirection.y=0;

	dragDirection.x=0;
	dragDirection.y=0;

	breakTime=0;
	lastTickCount=0;
}

void CycloDetector::startPan(int64 tickCount)
{
	//Start the pan
	dragStart=true;

	//Initialize a new stroke
	strokes.clear();
	stroke=new Stroke();
	stroke->setStartTick(tickCount);
	
	testPoint.x=touch.getResolutionX()/2;
	testPoint.y=touch.getResolutionY()/2;

	//Save the last tick
	lastTickCount=tickCount;

	//Save the first drag point
	stroke->addStrokePoint(touchPoints[0]);
	startDragPoint=touchPoints[0];
}

void CycloDetector::cycloZoom()
{
	Mat zoomImage=Mat::zeros(480,640,CV_8UC3);

	if(touchPoints.size()!=1)
	{
		cancelPan();
	}

	if(touchPoints.size()==1)
	{
		if(dragStart==false)
		{
			dragStart=true;
			stroke=new Stroke();
			stroke->addStrokePoint(touchPoints[0]);
		}
		else
		{

			stroke->addStrokePoint(touchPoints[0]);
			std::cout<<"stroke count "<<stroke->getStrokePoints().size()<<std::endl;
			Point startPoint(stroke->getStrokePoints()[0].x*640/touch.getResolutionX(),stroke->getStrokePoints()[0].y*480/touch.getResolutionY());
			Point endPoint(touchPoints[0].x*640/touch.getResolutionX(),touchPoints[0].y*480/touch.getResolutionY());

		    Point testPoint=startPoint*0.5+endPoint*0.5;

			line(zoomImage,startPoint,endPoint,Scalar(0,255,255));

			circle(zoomImage,testPoint,3,Scalar(0,255,0));
			for(int i=0;i<stroke->getStrokePoints().size();++i)
			{

			}
		}

		imshow("zoom",zoomImage);
	}
	
}

void CycloDetector::fitStroke2Ellispe(Stroke stroke)
{
	Mat S=Mat::zeros(6,6,CV_32FC1);

	for(int i=0;i<stroke.getStrokePoints().size();++i)
	{
		float x=stroke.getStrokePoints()[i].x;
		float y=stroke.getStrokePoints()[i].y;

		Mat D=Mat::zeros(1,6,CV_32FC1);
		D.at<float>(0,0)=x*x;
		D.at<float>(0,1)=x*y;
		D.at<float>(0,2)=y*y;
		D.at<float>(0,3)=x;
		D.at<float>(0,4)=y;
		D.at<float>(0,5)=1;

		S+=D.t()*D;
	}

	Mat C=Mat::zeros(6,6,CV_32FC1);
	C.at<float>(0,2)=2.0;
	C.at<float>(1,1)=-1.0;
	C.at<float>(2,0)=2.0;

	Mat eigenMatrix=S.inv()*C;
	Mat eigenValues=Mat::zeros(6,1,CV_32FC1);
	Mat eigenVectors=Mat::zeros(6,6,CV_32FC1);

	cv::eigen(eigenMatrix,true,eigenValues,eigenVectors);

	int maxIndex=0;
	int maxEigenValue=eigenValues.at<float>(0,0);
	for(int i=1;i<6;++i)
	{
		if(eigenValues.at<float>(i,0)>maxEigenValue)
		{
			maxIndex=i;
			maxEigenValue=eigenValues.at<float>(i,0);
		}
	}

	Mat A=eigenVectors.col(maxIndex);

	float b=A.at<float>(1,0)/2;
	float c=A.at<float>(2,0);
	float d=A.at<float>(3,0)/2;
	float f=A.at<float>(4,0)/2;
	float g=A.at<float>(5,0);
	float a=A.at<float>(0,0);

	float num=b*b-a*c;
	float x0=(c*d-b*f)/num;
	float y0=(a*f-b*d)/num;

	Point ellipseCenter(x0,y0);

	float rotationAngle=0.5*std::atan(2*b/(a-c));

	float up=2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g);
	float down1=(b*b-a*c)*( (c-a)*sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a));
    float down2=(b*b-a*c)*( (a-c)*sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a));

	float res1=sqrt(up/down1);
	float res2=sqrt(up/down2);



}