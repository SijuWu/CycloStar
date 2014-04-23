#include "StdAfx.h"
#include "Touch.h"
#include <iostream>
#include <set>
#include <map>
#include <cassert>
#include <functional>

 

namespace PQ_SDK_MultiTouch
{
static std::vector<TouchPoint> touchPointList;
//The first means the state, the last means the angle
//0:RotateState, 1:X coordinate of the anchor point, 2:Y coordinate of the anchor point
//3:X coordinate of the other point, 4:Y coordinate of the other point, 5:Rotate angle
//RotateState:0 nonValid, 1 Start, 2 Anticlock,3 Clock,4 End
static double rotateParameters[6];
//0: SplitState, 1:X coordinate of the first point, 2:Y coordinate of the first point
//3:X coordinate of the second point, 4:Y coordinate of the second point, 5:Finger distance, 6:Distance ratio
static double splitParameters[7];
static int resolutionX;
static int resolutionY;
//If touch input exists,this value is true
static bool touchExist;
	Touch::Touch(void)
	{
		memset(m_pf_on_tges,0, sizeof(m_pf_on_tges));
	}

	Touch::~Touch(void)
	{
		DisconnectServer();
	}

	int Touch::Init()
{
	int err_code = PQMTE_SUCCESS;
	
	// initialize the handle functions of gestures;
	InitFuncOnTG();
	// set the functions on server callback
	SetFuncsOnReceiveProc();
	// connect server
	cout << " connect to server..." << endl;
	if((err_code = ConnectServer()) != PQMTE_SUCCESS){
		cout << " connect server fail, socket error code:" << err_code << endl;
		return err_code;
	}
	// send request to server
	cout << " connect success, send request." << endl;
	TouchClientRequest tcq = {0};
	tcq.type = RQST_RAWDATA_ALL | RQST_GESTURE_ALL;
	if((err_code = SendRequest(tcq)) != PQMTE_SUCCESS){
		cout << " send request fail, error code:" << err_code << endl;
		return err_code;
	}
	////////////you can set the move_threshold when the tcq.type is RQST_RAWDATA_INSIDE;
	////send threshold
	//int move_threshold = 0;// 0 pixel, receuve all the touch points that touching in the windows area of this client;
	//if((err_code = SendThreshold(move_threshold)) != PQMTE_SUCCESS){
	//	cout << " send threadhold fail, error code:" << err_code << endl;
	//	return err_code;
	//}
	
	//////// you can set the resolution of the touch point(raw data) here;
	//// setrawdata_resolution
	//int maxX = 32768, maxY = 32768;
	//if((err_code= SetRawDataResolution(maxX, maxY)) != PQMTE_SUCCESS){
	//	cout << " set raw data resolution fail, error code:" << err_code << endl;
	//}
	////////////////////////
	//get server resolution
	if((err_code = GetServerResolution(OnGetServerResolution, NULL)) != PQMTE_SUCCESS){
		cout << " get server resolution fail,error code:" << err_code << endl;
		return err_code;
	}
	//
	// start receiving
	cout << " send request success, start recv." << endl;
	return err_code;
}

void Touch:: InitFuncOnTG()
{
	// initialize the call back functions of toucha gestures;
	m_pf_on_tges[TG_TOUCH_START] = &Touch::OnTG_TouchStart;
	m_pf_on_tges[TG_DOWN] = &Touch::OnTG_Down;
	m_pf_on_tges[TG_MOVE] = &Touch::OnTG_Move;
	m_pf_on_tges[TG_UP] = &Touch::OnTG_Up;

	m_pf_on_tges[TG_SECOND_DOWN] = &Touch::OnTG_SecondDown;
	m_pf_on_tges[TG_SECOND_UP] = &Touch::OnTG_SecondUp;

	m_pf_on_tges[TG_SPLIT_START] = &Touch::OnTG_SplitStart;
	m_pf_on_tges[TG_SPLIT_APART] = &Touch::OnTG_SplitApart;
	m_pf_on_tges[TG_SPLIT_CLOSE] = &Touch::OnTG_SplitClose;
	m_pf_on_tges[TG_SPLIT_END] = &Touch::OnTG_SplitEnd;

	m_pf_on_tges[TG_ROTATE_START] = &Touch::OnTG_RotateStart;
	m_pf_on_tges[TG_ROTATE_ANTICLOCK] = &Touch::OnTG_RotateAnticlock;
	m_pf_on_tges[TG_ROTATE_CLOCK] = &Touch::OnTG_RotateClock;
	m_pf_on_tges[TG_ROTATE_END] = &Touch::OnTG_RotateEnd;

	m_pf_on_tges[TG_TOUCH_END] = &Touch::OnTG_TouchEnd;

}
void Touch::SetFuncsOnReceiveProc()
{
	PFuncOnReceivePointFrame old_rf_func = SetOnReceivePointFrame(&Touch::OnReceivePointFrame,this);
	PFuncOnReceiveGesture old_rg_func = SetOnReceiveGesture(&Touch::OnReceiveGesture,this);
	PFuncOnServerBreak old_svr_break = SetOnServerBreak(&Touch::OnServerBreak,NULL);
	PFuncOnReceiveError old_rcv_err_func = SetOnReceiveError(&Touch::OnReceiveError,NULL);
	PFuncOnGetDeviceInfo old_gdi_func = SetOnGetDeviceInfo(&Touch::OnGetDeviceInfo,NULL);
}

void Touch:: OnReceivePointFrame(int frame_id, int time_stamp, int moving_point_count, const TouchPoint * moving_point_array, void * call_back_object)
{
	Touch * touch = static_cast<Touch*>(call_back_object);
	assert(touch != NULL);
	const char * tp_event[] = 
	{
		"down",
		"move",
		"up",
	};
	
	//cout << " frame_id:" << frame_id << " time:"  << time_stamp << " ms" << " moving point count:" << moving_point_count << endl;
	
	//cout<<moving_point_count<<endl;
	touchPointList.clear();
	for(int i = 0; i < moving_point_count; ++ i){
		TouchPoint tp = moving_point_array[i];
		touch->OnTouchPoint(tp);

		touchPointList.push_back(tp);
		if(i>=9)
			break;
	}

	

	//throw exception("test exception here");
}
void Touch:: OnReceiveGesture(const TouchGesture & ges, void * call_back_object)
{
	Touch * touch = static_cast<Touch*>(call_back_object);
	assert(touch != NULL);
	touch->OnTouchGesture(ges);
	//throw exception("test exception here");
}
void Touch:: OnServerBreak(void * param, void * call_back_object)
{
	// when the server break, disconenct server;
	//cout << "server break, disconnect here" << endl;
	DisconnectServer();
}
void Touch::OnReceiveError(int err_code, void * call_back_object)
{
	switch(err_code)
	{
	case PQMTE_RCV_INVALIDATE_DATA:
		cout << " error: receive invalidate data." << endl;
		break;
	case PQMTE_SERVER_VERSION_OLD:
		cout << " error: the multi-touch server is old for this client, please update the multi-touch server." << endl;
		break;
	case PQMTE_EXCEPTION_FROM_CALLBACKFUNCTION:
		cout << "**** some exceptions thrown from the call back functions." << endl;
		assert(0); //need to add try/catch in the callback functions to fix the bug;
		break;
	default:
		cout << " socket error, socket error code:" << err_code << endl;
	}
}
void Touch:: OnGetServerResolution(int x, int y, void * call_back_object)
{
	cout << " server resolution:" << x << "," << y << endl;
	resolutionX=x;
	resolutionY=y;
}
void Touch::OnGetDeviceInfo(const TouchDeviceInfo & deviceinfo,void *call_back_object)
{
	cout << " touch screen, SerialNumber: " << deviceinfo.serial_number <<",(" << deviceinfo.screen_width << "," << deviceinfo.screen_height << ")."<<  endl;
}
// here, just record the position of point,
//	you can do mouse map like "OnTG_Down" etc;
void Touch:: OnTouchPoint(const TouchPoint & tp)
{
	switch(tp.point_event)
	{
	case TP_DOWN:
	/*	cout << "  point " << tp.id << " come at (" << tp.x << "," << tp.y 
			<< ") width:" << tp.dx << " height:" << tp.dy << endl;*/
		break;
	case TP_MOVE:
		/*cout << "  point " << tp.id << " move at (" << tp.x << "," << tp.y 
			<< ") width:" << tp.dx << " height:" << tp.dy << endl;*/
		break;
	case TP_UP:
		/*cout << "  point " << tp.id << " leave at (" << tp.x << "," << tp.y 
			<< ") width:" << tp.dx << " height:" << tp.dy << endl;*/
		break;
	}
}
void Touch:: OnTouchGesture(const TouchGesture & tg)
{
	if(TG_NO_ACTION == tg.type)
		return ;
	
	assert(tg.type <= TG_TOUCH_END);
	DefaultOnTG(tg,this);
	PFuncOnTouchGesture pf = m_pf_on_tges[tg.type];
	if(NULL != pf){
		pf(tg,this);
	}
}
void Touch:: OnTG_TouchStart(const TouchGesture & tg,void * call_object)
{
	touchExist=true;
	assert(tg.type == TG_TOUCH_START);
	//cout << "  here, the touch start, initialize something." << endl;
}
void Touch:: DefaultOnTG(const TouchGesture & tg,void * call_object) // just show the gesture
{
	//cout <<"ges,name:"<< GetGestureName(tg) << " type:" << tg.type << ",param size:" << tg.param_size << " ";
	//for(int i = 0; i < tg.param_size; ++ i)
	//	cout << tg.params[i] << " ";
	//cout << endl;
}
void Touch:: OnTG_Down(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_DOWN && tg.param_size >= 2);
	//cout << "  the single finger touching at :( " 
		//<< tg.params[0] << "," << tg.params[1] << " )" << endl;
}
void Touch:: OnTG_Move(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_MOVE && tg.param_size >= 2);
	/*cout << "  the single finger moving on the screen at :( " 
		<< tg.params[0] << "," << tg.params[1] << " )" << endl;*/
}
void Touch:: OnTG_Up(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_UP && tg.param_size >= 2);
	/*cout << " the single finger is leaving the screen at :( " 
		<< tg.params[0] << "," << tg.params[1] << " )" << endl;*/
}
//
void Touch:: OnTG_SecondDown(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_SECOND_DOWN && tg.param_size >= 4);
	/*cout << "  the second finger touching at :( " 
		<< tg.params[0] << "," << tg.params[1] << " ),"
		<< " after the first finger touched at :( "
		<< tg.params[2] << "," << tg.params[3] << " )" << endl;*/
}
void Touch:: OnTG_SecondUp(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_SECOND_UP && tg.param_size >= 4);
	/*cout << "  the second finger is leaving at :( " 
		<< tg.params[0] << "," << tg.params[1] << " ),"
		<< " while the first finger still anchored around :( "
		<< tg.params[2] << "," << tg.params[3] << " )" << endl;*/
}
//
void Touch:: OnTG_SplitStart(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_SPLIT_START && tg.param_size >= 4);
	//cout << "  the two fingers is splitting with one finger at: ( " 
	//	<< tg.params[0] << "," << tg.params[1] << " ),"
	//	<< " , the other at :( "
	//	<< tg.params[2] << "," << tg.params[3] << " )" << endl;

	splitParameters[0]=1;
	splitParameters[1]=tg.params[0];
	splitParameters[2]=tg.params[1];
	splitParameters[3]=tg.params[2];
	splitParameters[4]=tg.params[3];
	splitParameters[5]=0;
}

void Touch:: OnTG_SplitApart(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_SPLIT_APART && tg.param_size >= 1);
	//cout << "  the two fingers is splitting apart with there distance incresed by " 
	//	<< tg.params[0]
	//	<< " with a ratio :" << tg.params[1]
	//	<< endl;

		splitParameters[0]=2;
		splitParameters[1]=tg.params[2];
		splitParameters[2]=tg.params[3];
		splitParameters[3]=tg.params[4];
		splitParameters[4]=tg.params[5];
		splitParameters[5]=tg.params[0];
		splitParameters[6]=tg.params[1];
}
void Touch:: OnTG_SplitClose(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_SPLIT_CLOSE && tg.param_size >= 1);
	//cout << "  the two fingers is splitting close with there distance decresed by " 
	//	<< tg.params[0]
	//	<< " with a ratio :" << tg.params[1]
	//	<< endl;

		splitParameters[0]=3;
		splitParameters[1]=tg.params[2];
		splitParameters[2]=tg.params[3];
		splitParameters[3]=tg.params[4];
		splitParameters[4]=tg.params[5];
		splitParameters[5]=tg.params[0];
		splitParameters[6]=tg.params[1];
}
void Touch:: OnTG_SplitEnd(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_SPLIT_END);
	//cout << "  the two splitting fingers with one finger at: ( " 
	//	<< tg.params[0] << "," << tg.params[1] << " ),"
	//	<< " , the other at :( "
	//	<< tg.params[2] << "," << tg.params[3] << " )" 
	//	<< " will end" << endl;

	splitParameters[0]=4;
		splitParameters[1]=tg.params[0];
		splitParameters[2]=tg.params[1];
		splitParameters[3]=tg.params[2];
		splitParameters[4]=tg.params[3];
		
}

void Touch::OnTG_RotateStart(const TouchGesture & tg,void * call_object)
{
	assert(tg.type==TG_ROTATE_START&& tg.param_size >= 4);
	rotateParameters[0]=1;
	rotateParameters[1]=tg.params[0];
	rotateParameters[2]=tg.params[1];
	rotateParameters[3]=tg.params[2];
	rotateParameters[4]=tg.params[3];
	rotateParameters[5]=0;
	//cout<<"Start"<<endl;
}
void Touch::OnTG_RotateAnticlock(const TouchGesture & tg,void * call_object)
{
	assert(tg.type==TG_ROTATE_ANTICLOCK	&& tg.param_size >= 5);
	rotateParameters[0]=2;
	rotateParameters[1]=tg.params[1];
	rotateParameters[2]=tg.params[2];
	rotateParameters[3]=tg.params[3];
	rotateParameters[4]=tg.params[4];
	rotateParameters[5]=tg.params[0];
	//cout<<"Anticlock"<<endl;
}
void Touch::OnTG_RotateClock(const TouchGesture & tg,void * call_object)
{
	assert(tg.type==TG_ROTATE_CLOCK&& tg.param_size >= 5);
	rotateParameters[0]=3;
	rotateParameters[1]=tg.params[1];
	rotateParameters[2]=tg.params[2];
	rotateParameters[3]=tg.params[3];
	rotateParameters[4]=tg.params[4];
	rotateParameters[5]=tg.params[0];
	//cout<<"Clock"<<endl;
}
void Touch::OnTG_RotateEnd(const TouchGesture & tg,void * call_object)
{
	assert(tg.type==TG_ROTATE_END&& tg.param_size >= 4);
	rotateParameters[0]=4;
	rotateParameters[1]=tg.params[0];
	rotateParameters[2]=tg.params[1];
	rotateParameters[3]=tg.params[2];
	rotateParameters[4]=tg.params[3];
	rotateParameters[5]=0;
	//cout<<"End"<<endl;
}

// OnTG_TouchEnd: to clear what need to clear
void Touch:: OnTG_TouchEnd(const TouchGesture & tg,void * call_object)
{
	assert(tg.type == TG_TOUCH_END);
	touchExist=false;
	//cout << "  all the fingers is leaving and there is no fingers on the screen." << endl;
}


/////////////////////////// functions ///////////////////////////////////
}

 std::vector<TouchPoint>& Touch::getTouchPointList()
{
	return touchPointList;
}

 int Touch::getResolutionX()
 {
	 return resolutionX;
 }

 int Touch::getResolutionY()
 {
	 return resolutionY;
 }

 double* Touch::getRotateParameters()
 {
	 return rotateParameters;
 }

 double* Touch::getSplitParameters()
 {
	 return splitParameters;
 }

 bool Touch::checkTouch()
 {
	 return touchExist;
 }