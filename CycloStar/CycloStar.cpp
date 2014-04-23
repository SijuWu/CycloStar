// CycloStar.cpp : définit le point d'entrée pour l'application console.
//

#include "stdafx.h"
#include "opencv2/opencv.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include "CycloDetector.h"
using namespace cv;
using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{
	CycloDetector cycloDetector;
	cycloDetector.init();
	while(true)
	{
		cycloDetector.detectTouch();
		//cycloDetector.cycloPan();
		cycloDetector.cycloZoom();

		Mat ellipse=Mat::zeros(480,640,CV_8UC3);
		vector<Point> ellipsePoints;
		vector<float>X;
		vector<float>Y;

		float mean_x=200;
		float mean_y=200;

		for(int i=0;i<360;++i)
		{
			ellipsePoints.push_back(Point(100*cos(CV_PI/180*i)+200,50*sin(CV_PI/180*i)+200));
			cv::circle(ellipse,Point(100*cos(CV_PI/180*i)+200,50*sin(CV_PI/180*i)+200),1,Scalar(0,255,0));
			X.push_back(100*cos(CV_PI/180*i));
			Y.push_back(50*sin(CV_PI/180*i));
		}

		Mat D=Mat::zeros(360,5,CV_32FC1);
		float* pointerD=(float*)D.data;

		for(int i=0;i<360;++i)
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
			for(int i=0;i<360;++i)
			{
				sumD.at<float>(0,j)=sumD.at<float>(0,j)+D.at<float>(i,j);
				
			}
			/*std::cout<<sumD.at<float>(1,j)<<std::endl;*/
		}
	
		
		Mat A=sumD*(D.t()*D).inv();
		std::cout<<A<<std::endl;

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

			float mean_x=cos_phi*mean_xtemp - sin_phi*mean_ytemp;
			float mean_y=sin_phi*mean_xtemp + cos_phi*mean_ytemp;
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
				std::cout<<X0<<" "<<Y0<<std::endl;
				float F=1 + (d*d)/(4*a) + (e*e)/(4*c);
				a=sqrt(F/a);
				b=sqrt(F/c);

				float long_axis=2*max(a,b);
				float short_axis=2*min(a,b);

				std::cout<<long_axis<<std::endl;
				std::cout<<short_axis<<std::endl;
			
		}
		//Mat S=Mat::zeros(6,6,CV_32FC1);
		//Mat D=Mat::zeros(360,6,CV_32FC1);
		//float* pointerD=(float*)D.data;

		//for(int i=0;i<360;++i)
		//{
		//	float x=100*cos(CV_PI/180*i)+200;
		//	float y=50*sin(CV_PI/180*i)+200;

		//	D.at<float>(i,0)=x*x;
		//	D.at<float>(i,1)=x*y;
		//	D.at<float>(i,2)=y*y;
		//	D.at<float>(i,3)=x;
		//	D.at<float>(i,4)=y;
		//	D.at<float>(i,5)=1;
		//	//*(pointerD+i*6)=x*x;
		//	//*(pointerD+i*6+1)=x*y;
		//	//*(pointerD+i*6+2)=y*y;
		//	//*(pointerD+i*6+3)=x;
		//	//*(pointerD+i*6+4)=y;
		//	//*(pointerD+i*6+5)=1;
		//}

		//S=D.t()*D;
		///*for(int i=0;i<360;++i)
		//{
		//	float x=100*cos(CV_PI/180*i)+200;
		//	float y=50*sin(CV_PI/180*i)+200;

		//	Mat D=Mat::zeros(1,6,CV_32FC1);
		//	D.at<float>(0,0)=x*x;
		//	D.at<float>(0,1)=x*y;
		//	D.at<float>(0,2)=y*y;
		//	D.at<float>(0,3)=x;
		//	D.at<float>(0,4)=y;
		//	D.at<float>(0,5)=1;

		//	S+=D.t()*D;
		//}*/


		//Mat C=Mat::zeros(6,6,CV_32FC1);
		//C.at<float>(0,2)=2.0;
		//C.at<float>(1,1)=-1.0;
		//C.at<float>(2,0)=2.0;


		//
		//Mat eigenMatrix=S.inv()*C;


		//Mat eigenValues=Mat::zeros(6,1,CV_32FC1);
		//Mat eigenVectors=Mat::zeros(6,6,CV_32FC1);

		//cv::eigen(eigenMatrix,true,eigenValues,eigenVectors);

		//int maxIndex=0;
		//int maxEigenValue=eigenValues.at<float>(0,0);
		//for(int i=1;i<6;++i)
		//{
		//	if(eigenValues.at<float>(i,0)>maxEigenValue)
		//	{
		//		maxIndex=i;
		//		maxEigenValue=eigenValues.at<float>(i,0);
		//	}
		//}
		//std::cout<<"eigen value "<<eigenValues<<std::endl;
		//Mat A=eigenVectors.col(maxIndex);

		//float b=A.at<float>(1,0)/2;
		//float c=A.at<float>(2,0);
		//float d=A.at<float>(3,0)/2;
		//float f=A.at<float>(4,0)/2;
		//float g=A.at<float>(5,0);
		//float a=A.at<float>(0,0);

		//float num=b*b-a*c;
		//float x0=(c*d-b*f)/num;
		//float y0=(a*f-b*d)/num;

		//Point ellipseCenter(x0,y0);

		//float rotationAngle=0.5*std::atan(2*b/(a-c));

		//float up=2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g);
		//float down1=(b*b-a*c)*( (c-a)*sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a));
		//float down2=(b*b-a*c)*( (a-c)*sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a));

		//float res1=sqrt(up/down1);
		//float res2=sqrt(up/down2);

		//std::cout<<"center "<<ellipseCenter<<std::endl;

		imshow("ellipse",ellipse);

		if(waitKey(1)=='q')
		{
			break;
		}

	}




	return 0;
}

