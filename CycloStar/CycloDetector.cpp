#include "StdAfx.h"
#include "CycloDetector.h"

static float scale=1;
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

vector<TouchPoint>& Stroke::getStrokePoints()
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

void Stroke::setEllipseParameters(float X0,float Y0, float a,float b,float long_axis, float short_axis,float cos_phi,float sin_phi)
{
	this->X0=X0;
	this->Y0=Y0;

	this->a=a;
	this->b=b;

	this->long_axis=long_axis;
	this->short_axis=short_axis;

	this->cos_phi=cos_phi;
	this->sin_phi=sin_phi;
}

float Stroke::getX0()
{
	return this->X0;
}

float Stroke::getY0()
{
	return this->Y0;
}

float Stroke::getA()
{
	return this->a;
}

float Stroke::getB()
{
	return this->b;
}

float Stroke::getCos_phi()
{
	return this->cos_phi;
}

float Stroke::getSin_phi()
{
	return this->sin_phi;
}
float Stroke::getLong_axis()
{
	return this->long_axis;
}

float Stroke::getShort_axis()
{
	return this->short_axis;
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
	stroke=new Stroke();
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
		stroke->getStrokePoints().clear();
	}

	if(touchPoints.size()==1)
	{
		stroke->addStrokePoint(touchPoints[0]);


		for(int i=0;i<stroke->getStrokePoints().size();++i)
		{
			Point point(stroke->getStrokePoints()[i].x*640/touch.getResolutionX(),stroke->getStrokePoints()[i].y*480/touch.getResolutionY());
			circle(zoomImage,point,1,Scalar(0,255,0));
		}	
	}
	imshow("zoom",zoomImage);
	
	fitStroke2Ellispe(*stroke);
}

void CycloDetector::cycloCheck()
{
	Mat strokeEllipse=Mat::zeros(480,640,CV_8UC3);
	if(touchPoints.size()!=1)
	{
		stroke->getStrokePoints().clear();
		strokes.clear();
		return;
	}
	else
	{
		stroke->getStrokePoints().push_back(touchPoints[0]);

		int pointCount=stroke->getStrokePoints().size();

		for(int i=0;i<stroke->getStrokePoints().size();++i)
		{
			circle(strokeEllipse,Point(stroke->getStrokePoints()[i].x*640/touch.getResolutionX(),stroke->getStrokePoints()[i].y*480/touch.getResolutionY()),1,Scalar(0,255,255));
		}

		if(strokes.size()!=0)
		{
			for(int i=0;i<strokes[strokes.size()-1]->getStrokePoints().size();++i)
			{
				circle(strokeEllipse,Point(strokes[strokes.size()-1]->getStrokePoints()[i].x*640/touch.getResolutionX(),strokes[strokes.size()-1]->getStrokePoints()[i].y*480/touch.getResolutionY()),1,Scalar(0,255,255));
			}
		}

		if(strokes.size()==0&&stroke->getStrokePoints().size()<5)
			return;

		TouchPoint startPoint=stroke->getStrokePoints()[0];
		TouchPoint lastEndPoint=stroke->getStrokePoints()[pointCount-2];
		TouchPoint endPoint=stroke->getStrokePoints()[pointCount-1];

		float lastDis=getDistance(startPoint,lastEndPoint);
		float dis=getDistance(startPoint,endPoint);

		//When the stroke reverse, fit the stroke to ellipse
		if((lastDis-dis)>1)
		{
			strokes.push_back(stroke);
			stroke=new Stroke();
			stroke->addStrokePoint(touchPoints[0]);

			if(strokes.size()>=2)
			{
				Stroke* stroke1=strokes[strokes.size()-2];
				Stroke* stroke2=strokes[strokes.size()-1];

				Stroke* stroke3=new Stroke();
				stroke3->getStrokePoints().insert(stroke3->getStrokePoints().end(),stroke1->getStrokePoints().begin(),stroke1->getStrokePoints().end());
				stroke3->getStrokePoints().insert(stroke3->getStrokePoints().end(),stroke2->getStrokePoints().begin(),stroke2->getStrokePoints().end());

				fitStroke2Ellispe(*stroke3,*stroke2);
				fitStroke2Circle(*stroke3,*stroke2);
				float powA=powf(stroke2->getLong_axis()/2,2);
				float powB=powf(stroke2->getShort_axis()/2,2);

				float eccentricity=sqrtf(1.0f-powB/powA);
				std::cout<<"eccentricity "<<eccentricity<<" a "<<stroke2->getA()<<" b "<<stroke2->getB()<<" long "<<stroke2->getLong_axis()<<" short "<<stroke2->getShort_axis()<<std::endl;
			}			
		}

		if(strokes.size()>=2)
		{
			Stroke* lastStroke=strokes[strokes.size()-1];

			cv::Matx22f R(lastStroke->getCos_phi(),lastStroke->getSin_phi(),
				         -lastStroke->getSin_phi(),lastStroke->getCos_phi());

			for(int i=0;i<360;++i)
			{
				float xinit=lastStroke->getX0()+ lastStroke->getA()*cos( CV_PI/180*i );
				float yinit=lastStroke->getY0()+ lastStroke->getB()*sin( CV_PI/180*i );

				cv::Matx21f initialPoint(xinit,yinit);
				cv::Matx21f rotatedPoint=R*initialPoint;
				cv::circle(strokeEllipse,Point(rotatedPoint(0)*640/touch.getResolutionX(),rotatedPoint(1)*480/touch.getResolutionY()),1,Scalar(0,0,255));
				
			}

			cv::Matx21f initialCenter(lastStroke->getX0(),lastStroke->getY0());
			cv::Matx21f rotatedCenter=R*initialCenter;
			cv::circle(strokeEllipse,Point(rotatedCenter(0)*640/touch.getResolutionX(),rotatedCenter(1)*480/touch.getResolutionY()),1,Scalar(0,0,255));

			TouchPoint startPoint1=stroke->getStrokePoints()[0];
			int curveCenterIndex=floorf(stroke->getStrokePoints().size()/2);
			TouchPoint startPoint2=stroke->getStrokePoints()[curveCenterIndex];

			Vec2f lastVec(startPoint1.x-rotatedCenter(0),startPoint1.y-rotatedCenter(1));
			Vec2f vec(startPoint2.x-rotatedCenter(0),startPoint2.y-rotatedCenter(1));

			float cross=lastVec(0)*vec(1)-lastVec(1)*vec(0);

			if(stroke->getStrokePoints().size()>=2)
			{
				TouchPoint lastPoint=stroke->getStrokePoints()[stroke->getStrokePoints().size()-2];
				TouchPoint point=stroke->getStrokePoints()[stroke->getStrokePoints().size()-1];

				float lastDis1=(lastPoint.x-rotatedCenter(0))*(lastPoint.x-rotatedCenter(0))+(lastPoint.y-rotatedCenter(1))*(lastPoint.y-rotatedCenter(1));
				float dis1=(point.x-rotatedCenter(0))*(point.x-rotatedCenter(0))+(point.y-rotatedCenter(1))*(point.y-rotatedCenter(1));
				float motion=(lastPoint.x-point.x)*(lastPoint.x-point.x)+(lastPoint.y-point.y)*(lastPoint.y-point.y);
				lastDis1=sqrt(lastDis1);
				dis1=sqrt(dis1);
				motion=sqrt(motion);

				float averageDis=(lastDis1+dis1)/2;

				float angleChange=motion/averageDis;

				std::cout<<"anglevariation "<<angleChange<<std::endl;
				float K;
				if(cross>0)
				  K=0.14;
				if(cross<0)
					K=-0.32;
				std::cout<<"range"<<1+K*angleChange<<std::endl;

				scale=scale*(1+K*angleChange);
				std::cout<<"scale "<<scale<<std::endl;

			}
			//if((lastDis-dis)<0)
			//{
			//	if(cross>0)
			//		std::cout<<"clockwise "<<cross<<std::endl;
			//	if(cross<0)
			//		std::cout<<"anticlockwise "<<cross<<std::endl;
			//}
			
		}
	}

	imshow("strokeEllipse",strokeEllipse);
}

void CycloDetector::fitStroke2Ellispe(Stroke& stroke)
{
	Mat resultImage=Mat::zeros(480,640,CV_8UC3);
	vector<float>X;
	vector<float>Y;

	if(stroke.getStrokePoints().size()<5)
		return;
	for(int i=0;i<stroke.getStrokePoints().size();++i)
	{
		float x=stroke.getStrokePoints()[i].x;
		float y=stroke.getStrokePoints()[i].y;

		X.push_back(x);
		Y.push_back(y);
		Mat D=Mat::zeros(1,6,CV_32FC1);
	}

	float mean_x=mean(X);
	float mean_y=mean(Y);

	reduce(X,mean_x);
	reduce(Y,mean_y);

	Mat D=Mat::zeros(stroke.getStrokePoints().size(),5,CV_32FC1);
	float* pointerD=(float*)D.data;

	for(int i=0;i<stroke.getStrokePoints().size();++i)
	{
		D.at<float>(i,0)=X[i]*X[i];
		D.at<float>(i,1)=X[i]*Y[i];
		D.at<float>(i,2)=Y[i]*Y[i];
		D.at<float>(i,3)=X[i];
		D.at<float>(i,4)=Y[i];
	}

	Mat sumD=Mat::zeros(1,5,CV_32FC1);
	for(int j=0;j<5;++j)
	{
		for(int i=0;i<stroke.getStrokePoints().size();++i)
		{
			sumD.at<float>(0,j)=sumD.at<float>(0,j)+D.at<float>(i,j);		
		}
	}


	Mat A=sumD*(D.t()*D).inv();

	float a=A.at<float>(0,0);
	float b=A.at<float>(0,1);
	float c=A.at<float>(0,2);
	float d=A.at<float>(0,3);
	float e=A.at<float>(0,4);

	float orientation_tolerance=0.001;

	float orientation_rad;
	float cos_phi;
	float sin_phi;

	if(std::min(abs(b/a),abs(b/c))>orientation_tolerance)
		{
			orientation_rad=0.5*std::atan(b/(c-a));
			cos_phi=cos(orientation_rad);
			sin_phi=sin(orientation_rad);
			
			float atemp=a;
			float btemp=b;
			float ctemp=c;
			float dtemp=d;
			float etemp=e;
			a= atemp*cos_phi*cos_phi - btemp*cos_phi*sin_phi + ctemp*sin_phi*sin_phi;
			b=0;
			c= atemp*sin_phi*sin_phi + btemp*cos_phi*sin_phi + ctemp*cos_phi*cos_phi;
			d= dtemp*cos_phi - etemp*sin_phi;
			e= dtemp*sin_phi + etemp*cos_phi;

			float mean_xtemp=mean_x;
			float mean_ytemp=mean_y;

			mean_x=cos_phi*mean_xtemp - sin_phi*mean_ytemp;
			mean_y=sin_phi*mean_xtemp + cos_phi*mean_ytemp;
		}
		else
		{
			orientation_rad=0;
			cos_phi= cos( orientation_rad );
			sin_phi=sin(orientation_rad);
		}

		float test=a*c;

		if(test>0)
		{
			if(a<0)
			{
				a=-a;
				c=-c;
				d=-d;
				e=-e;
				
			}

				float X0=mean_x-d/2/a;
				float Y0=mean_y-e/2/c;
				
				float F=1 + (d*d)/(4*a) + (e*e)/(4*c);
				a=sqrt(F/a);
				b=sqrt(F/c);

				float long_axis=2*max(a,b);
				float short_axis=2*min(a,b);
				//std::cout<<"X0 "<<X0<<" Y0 "<<Y0<<" a "<<a<<" b "<<b<<" long "<<long_axis<<" short "<<short_axis<<" angle "<<orientation_rad<<std::endl;				

				stroke.setEllipseParameters(X0,Y0,a,b,long_axis,short_axis,cos_phi,sin_phi);
				cv::Matx22f R(cos_phi,sin_phi,-sin_phi,cos_phi);
			
				for(int i=0;i<360;++i)
				{
					float xinit=X0 + a*cos( CV_PI/180*i );
					float yinit=Y0 + b*sin( CV_PI/180*i );

					cv::Matx21f initialPoint(xinit,yinit);
					cv::Matx21f rotatedPoint=R*initialPoint;

					cv::circle(resultImage,Point(rotatedPoint(0)*640/touch.getResolutionX(),rotatedPoint(1)*480/touch.getResolutionY()),1,Scalar(0,0,255));
				}
		}

		else
			std::cout<<"no ellipse fitting"<<std::endl;
		imshow("result",resultImage);
}

void CycloDetector::fitStroke2Ellispe(Stroke& stroke1,Stroke& stroke2)
{
	Mat resultImage=Mat::zeros(480,640,CV_8UC3);
	vector<float>X;
	vector<float>Y;

	if(stroke1.getStrokePoints().size()<5)
		return;
	for(int i=0;i<stroke1.getStrokePoints().size();++i)
	{
		float x=stroke1.getStrokePoints()[i].x;
		float y=stroke1.getStrokePoints()[i].y;

		X.push_back(x);
		Y.push_back(y);
		Mat D=Mat::zeros(1,6,CV_32FC1);
	}

	float mean_x=mean(X);
	float mean_y=mean(Y);

	reduce(X,mean_x);
	reduce(Y,mean_y);

	Mat D=Mat::zeros(360,5,CV_32FC1);
	float* pointerD=(float*)D.data;

	for(int i=0;i<stroke1.getStrokePoints().size();++i)
	{
		D.at<float>(i,0)=X[i]*X[i];
		D.at<float>(i,1)=X[i]*Y[i];
		D.at<float>(i,2)=Y[i]*Y[i];
		D.at<float>(i,3)=X[i];
		D.at<float>(i,4)=Y[i];
	}

	Mat sumD=Mat::zeros(1,5,CV_32FC1);
	for(int j=0;j<5;++j)
	{
		for(int i=0;i<stroke1.getStrokePoints().size();++i)
		{
			sumD.at<float>(0,j)=sumD.at<float>(0,j)+D.at<float>(i,j);		
		}
	}


	Mat A=sumD*(D.t()*D).inv();

	float a=A.at<float>(0,0);
	float b=A.at<float>(0,1);
	float c=A.at<float>(0,2);
	float d=A.at<float>(0,3);
	float e=A.at<float>(0,4);

	float orientation_tolerance=0.001;

	float orientation_rad;
	float cos_phi;
	float sin_phi;

	if(std::min(abs(b/a),abs(b/c))>orientation_tolerance)
		{
			orientation_rad=0.5*std::atan(b/(c-a));
			cos_phi=cos(orientation_rad);
			sin_phi=sin(orientation_rad);
			
			float atemp=a;
			float btemp=b;
			float ctemp=c;
			float dtemp=d;
			float etemp=e;
			a= atemp*cos_phi*cos_phi - btemp*cos_phi*sin_phi + ctemp*sin_phi*sin_phi;
			b=0;
			c= atemp*sin_phi*sin_phi + btemp*cos_phi*sin_phi + ctemp*cos_phi*cos_phi;
			d= dtemp*cos_phi - etemp*sin_phi;
			e= dtemp*sin_phi + etemp*cos_phi;

			float mean_xtemp=mean_x;
			float mean_ytemp=mean_y;

			mean_x=cos_phi*mean_xtemp - sin_phi*mean_ytemp;
			mean_y=sin_phi*mean_xtemp + cos_phi*mean_ytemp;
		}
		else
		{
			orientation_rad=0;
			cos_phi= cos( orientation_rad );
			sin_phi=sin(orientation_rad);
		}

		float test=a*c;

		if(test>0)
		{
			if(a<0)
			{
				a=-a;
				c=-c;
				d=-d;
				e=-e;
				
			}

				float X0=mean_x-d/2/a;
				float Y0=mean_y-e/2/c;
				
				float F=1 + (d*d)/(4*a) + (e*e)/(4*c);
				a=sqrt(F/a);
				b=sqrt(F/c);

				float long_axis=2*max(a,b);
				float short_axis=2*min(a,b);
				//std::cout<<"X0 "<<X0<<" Y0 "<<Y0<<" a "<<a<<" b "<<b<<" long "<<long_axis<<" short "<<short_axis<<" angle "<<orientation_rad<<std::endl;				

				stroke2.setEllipseParameters(X0,Y0,a,b,long_axis,short_axis,cos_phi,sin_phi);
				cv::Matx22f R(cos_phi,sin_phi,-sin_phi,cos_phi);
			
				for(int i=0;i<360;++i)
				{
					float xinit=X0 + a*cos( CV_PI/180*i );
					float yinit=Y0 + b*sin( CV_PI/180*i );

					cv::Matx21f initialPoint(xinit,yinit);
					cv::Matx21f rotatedPoint=R*initialPoint;

					cv::circle(resultImage,Point(rotatedPoint(0)*640/touch.getResolutionX(),rotatedPoint(1)*480/touch.getResolutionY()),1,Scalar(0,0,255));
				}
		}

		else
			std::cout<<"no ellipse fitting"<<std::endl;
		imshow("result",resultImage);
}

void CycloDetector::fitStroke2Circle(Stroke& stroke1,Stroke& stroke2)
{
	Mat circleImage=Mat::zeros(480,640,CV_8UC3);
	vector<float>X;
	vector<float>Y;

	if(stroke1.getStrokePoints().size()<5)
		return;
	for(int i=0;i<stroke1.getStrokePoints().size();++i)
	{
		float x=stroke1.getStrokePoints()[i].x;
		float y=stroke1.getStrokePoints()[i].y;

		X.push_back(x);
		Y.push_back(y);
		Mat D=Mat::zeros(1,6,CV_32FC1);
	}

	float mean_x=mean(X);
	float mean_y=mean(Y);

	reduce(X,mean_x);
	reduce(Y,mean_y);

	Mat D=Mat::zeros(360,3,CV_32FC1);
	float* pointerD=(float*)D.data;

	for(int i=0;i<stroke1.getStrokePoints().size();++i)
	{
		D.at<float>(i,0)=X[i]*X[i]+Y[i]*Y[i];
		D.at<float>(i,1)=X[i];
		D.at<float>(i,2)=Y[i];
	}

	Mat sumD=Mat::zeros(1,3,CV_32FC1);
	for(int j=0;j<4;++j)
	{
		for(int i=0;i<stroke1.getStrokePoints().size();++i)
		{
			sumD.at<float>(0,j)=sumD.at<float>(0,j)+D.at<float>(i,j);		
		}
	}


	Mat A=sumD*(D.t()*D).inv();
	std::cout<<"A "<<A<<std::endl;
	float a=A.at<float>(0,0);
	float b=A.at<float>(0,1);
	float c=A.at<float>(0,2);

	float radius=sqrtf(1/a);
	float X0=mean_x-b/(2*a);
	float Y0=mean_y-c/(2*a);

	circle(circleImage,Point(X0*640/touch.getResolutionX(),Y0*480/touch.getResolutionY()),radius,Scalar(0,255,0));
	imshow("circleImage",circleImage);
}

float CycloDetector::mean(vector<float> coordinates)
{
	float mean=0;
	for(int i=0;i<coordinates.size();++i)
	{
		mean+=coordinates[i];
	}

	return mean/=coordinates.size();
}

void CycloDetector::reduce(vector<float>& coordinates, float value)
{
	for(int i=0;i<coordinates.size();++i)
	{
		coordinates[i]-=value;
	}
}