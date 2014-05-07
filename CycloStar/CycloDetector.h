#pragma once
#include "opencv2/opencv.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include "Touch.h"
#include <vector>

using namespace cv;
using namespace std;

class Stroke
{
public:
	Stroke();
	~Stroke();

	vector<TouchPoint>& getStrokePoints();
	void addStrokePoint(TouchPoint point);

	Point2f getStrokeDirection();
	void setStrokeDirection(Point2f direction);

	Point2f getStrokeSpeed();
	void setStrokeSpeed(Point2f speed);

	void setStartTick(int64 startTick);
	void setEndTick(int64 endTick);

	void setFrequence(float frequence);
	float getFrequence();

	void setGain(vector<Stroke*> strokes);
	float getGain();

	void setEllipseParameters(float X0,float Y0, float a,float b,float long_axis, float short_axis,float cos_phi,float sin_phi);

	float getX0();
	float getY0();
	float getA();
	float getB();
	float getLong_axis();
	float getShort_axis();
	float getCos_phi();
	float getSin_phi();



private:
	vector<TouchPoint> strokePoints;
	Point2f direction;
	Point2f speed;
	int64 startTick;
	int64 endTick;
	float frequence;
	float gain;

	float a;
	float b;

	float long_axis;
	float short_axis;

	

	float X0;
	float Y0;

	float cos_phi;
	float sin_phi;

	

};

class CycloDetector
{
public:
	CycloDetector(void);
	~CycloDetector(void);
	//Detect touch input
	void detectTouch();
	//Initialize touch
	void init();
	//Get the distance between two points
	float getDistance(TouchPoint point1,TouchPoint point2);
	//Get the direction of the stroke
	Point2f getDirection(TouchPoint endPoint,TouchPoint startPoint);
	//Check reverse
	bool checkReverse(Point2f direction,Point2f lastDirection);

	//Cyclo pan
	void cycloPan();
	void cancelPan();
	//Cyclo zoom
	void cycloZoom();

	void cycloCheck();

	void startPan(int64 tickCount);

	void fitStroke2Ellispe(Stroke& stroke);
	void fitStroke2Ellispe(Stroke& stroke1,Stroke&stroke2);
	void fitStroke2Circle(Stroke& stroke1,Stroke&stroke2);
	//Get mean value of a vector
	float mean(vector<float> coordinates);
	//Reduce a value for each element 
	void reduce(vector<float>& coordinates, float value);
private:
	Touch touch;
	//The list of touch points of the last frame
	vector<TouchPoint> lastTouchPoints;
	vector<TouchPoint> touchPoints;

	TouchPoint startDragPoint;
	bool dragStart;
	bool firstDragEnd;
	Point2f lastDragDirection;
	Point2f dragDirection;
	Point2f firstDragDirection;
	float dragAmplitude;
	//Vector of strokes
	vector<Stroke*> strokes;
	//Active stroke
	Stroke* stroke;

	Point2f testPoint;

	//Break time of dragging
	float breakTime;
	int64 lastTickCount;
};

