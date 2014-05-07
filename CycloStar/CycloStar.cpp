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


float mean(vector<float> coordinates)
{
	float mean=0;
	for(int i=0;i<coordinates.size();++i)
	{
		mean+=coordinates[i];
	}

	return mean/=coordinates.size();
}

void reduce(vector<float>& coordinates, float value)
{
	for(int i=0;i<coordinates.size();++i)
	{
		coordinates[i]-=value;
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	CycloDetector cycloDetector;
	cycloDetector.init();
	while(true)
	{
		cycloDetector.detectTouch();
		//cycloDetector.cycloPan();
		//cycloDetector.cycloZoom();
		cycloDetector.cycloCheck();

		Mat ellipse=Mat::zeros(480,640,CV_8UC3);
		vector<Point> ellipsePoints;
		vector<Point> ellipseRPoints;
		vector<float>X;
		vector<float>Y;
		vector<float>Xr;
		vector<float>Yr;
		int m=180;

		float xcenter=200;
		float ycenter=200;

		float aInit=100;
		float bInit=50;
		float angle=45*CV_PI/180;

		for(int i=0;i<m;++i)
		{
			
			float x=aInit*cos(CV_PI/180*i)+xcenter;
			float y=bInit*sin(CV_PI/180*i)+ycenter;
			X.push_back(x);
			Y.push_back(y);

			float xr=(x-xcenter)*cos(angle)-(y-ycenter)*sin(angle)+xcenter;
			float yr=(x-xcenter)*sin(angle)+(y-ycenter)*cos(angle)+ycenter;
			Xr.push_back(xr);
			Yr.push_back(yr);

			ellipsePoints.push_back(Point(x,y));
			cv::circle(ellipse,Point(x,y),1,Scalar(0,255,0));
			ellipseRPoints.push_back(Point(xr,yr));
			cv::circle(ellipse,Point(xr,yr),1,Scalar(0,255,255));


		}

		float mean_x=mean(Xr);
		float mean_y=mean(Yr);

		reduce(Xr,mean_x);
		reduce(Yr,mean_y);

		Mat D=Mat::zeros(360,5,CV_32FC1);
		float* pointerD=(float*)D.data;

		for(int i=0;i<m;++i)
		{
			D.at<float>(i,0)=Xr[i]*Xr[i];
			D.at<float>(i,1)=Xr[i]*Yr[i];
			D.at<float>(i,2)=Yr[i]*Yr[i];
			D.at<float>(i,3)=Xr[i];
			D.at<float>(i,4)=Yr[i];
		}

		//std::cout<<"D'*D "<<D.t()*D<<std::endl;
		Mat sumD=Mat::zeros(1,5,CV_32FC1);
		for(int j=0;j<5;++j)
		{
			for(int i=0;i<m;++i)
			{
				sumD.at<float>(0,j)=sumD.at<float>(0,j)+D.at<float>(i,j);		
			}
		}
	
		
		Mat A=sumD*(D.t()*D).inv();
	
		//std::cout<<A<<std::endl;

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

				cv::Matx22f R(cos_phi,sin_phi,-sin_phi,cos_phi);
			
				for(int i=0;i<360;++i)
				{
					float xinit=X0 + a*cos( CV_PI/180*i );
					float yinit=Y0 + b*sin( CV_PI/180*i );

					cv::Matx21f initialPoint(xinit,yinit);
					cv::Matx21f rotatedPoint=R*initialPoint;

					cv::circle(ellipse,Point(rotatedPoint(0),rotatedPoint(1)),1,Scalar(0,0,255));
				}
		}
		

		imshow("ellipse",ellipse);

		if(waitKey(1)=='q')
		{
			break;
		}

	}




	return 0;
}

