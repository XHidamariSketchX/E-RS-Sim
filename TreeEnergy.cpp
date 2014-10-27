//============================================================================
// Name        : TreeEnergy.cpp
// Author      : XUJI
// Version     : 0.1
// Copyright   : 
// Description : 
//============================================================================

#include "TreeEnergy.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream.h >
#include <iomanip.h>

#define random(x) (rand()%x)

int choice=0;
int PK_Sim=50;
double PK_Cal=0.5;
Tree tree;

int root=0;
double energy=0;
double arrivalCount=0;

/********************Tools***********************************/
double mulp(double i){
	double m=1;
	int c=1;
	for(c=1;c<=i;c++){
		m*=c;
	}
	return m;
	cout<<"endl"<<endl;
}
double CC(double x,double y){
	double c=0;
	c=mulp(x)/(mulp(y)*mulp(x-y));

	return c;	
}

double Etx(double k,double d){
	if(k==0)return 0;
	return k*(50+0.01*d*d);
}
double Erx(double k){
	if(k==0)return 0;
	return k*50; 
}
/********************End of Tools***********************************/
/********************create tree**************************************/
double cal_distance(int x0,int y0,int x1,int y1){
	return sqrt(pow(x1-x0,2)+pow(y1-y0,2));
}
int check_conn(Tree *t,int x,int y){
	//@return father node index; 0 if Non of them is connected
	//@input x,y is the node location

	//visit all node s and find the closest node index then return it
	// check if  it is out of degree
	int size;
	int tmp_x,tmp_y;
	int cloest=-1;
	double distance;
	double dis_tmp=0;
	size=t->size;
	//initial distance to root;
	tmp_x=t->treeNode[t->root].x;
	tmp_y=t->treeNode[t->root].y;
	distance=cal_distance(tmp_x,tmp_y,x,y);

	if(distance<=DIS){
	//can connect to root;
		cloest=t->root;
	}
	//cout<<"dis=="<<distance<<endl;

	for(int i=1;i<size;i++){
		tmp_x=t->treeNode[i].x;
		tmp_y=t->treeNode[i].y;
		//calculate distance to next node
		dis_tmp=cal_distance(tmp_x,tmp_y,x,y);
		//cout<<"dis=="<<dis_tmp<<endl;
		//check if
		if(dis_tmp>DIS){
		// can not connect to current node
		}else{
		// it can connect to
			if(dis_tmp<distance){
			//current is closer
				distance=dis_tmp;
				cloest=i;
			}
		}
	}
	//cout<<"cloest=="<<cloest<<endl;
	return cloest;
}
void appendNode(Tree *t, int father,int current_node,int x,int y){
	CSNode *p;
	p=&t->treeNode[father];

	t->treeNode[current_node].father=p;
	if(p->firstchild==NULL){
	//if(t->treeNode[father].firstchild==0){
		//the node index has no child ,current_node is the first one
		p->firstchild=&t->treeNode[current_node];
		//cout<<"p->firstchild==NULL"<<endl;
	}else{
		//visit all the sibling from first one
		p=(p->firstchild);

		//cout<<"p== "<<p->nextsibling<<endl;

		while(p->nextsibling){
		//move to last sibling
			p=p->nextsibling;
		}
		//append to father
		p->nextsibling=&(t->treeNode[current_node]);
	}
}
void creatTree(Tree *t,int size){
	int x,y;
	int current_size=1;
	int father=0;
	if(size>MAX_TREE_SIZE)return;
	for(int i=1;i<size;i++){//append node 1-size to tree t
		//random a position
		x=random(LOC_SIZE);
		y=random(LOC_SIZE);

		while((father=check_conn(t,x,y))==-1){
			x=random(LOC_SIZE);
			y=random(LOC_SIZE);
		}
		//put node on the map;
		t->treeNode[i].x=x;
		t->treeNode[i].y=y;
		t->treeNode[i].dis_to_father=cal_distance(x,y,t->treeNode[father].x,t->treeNode[father].y);

		//set data
		t->treeNode[i].data=i;
		//end of set data
		//check connect again
		if(t->treeNode[i].dis_to_father<=DIS){
		//connected
			//cout<<"connected"<<endl;
		}else{
			cout<<"ERROR!! find wrong point"<<endl;
			cout<<"ERROR distance is "<<t->treeNode[i].dis_to_father<<endl;
		}
		//cout<<"node "<<i<<"'s father is "<<father<<endl;
		appendNode(t,father,i,x,y);
		t->size++;
	}
	if(size!=t->size){
		//error
	}
}
void initTree(Tree *t){
	t->root=0;
	t->treeNode[0].x=ROOT_LOC;
	t->treeNode[0].y=ROOT_LOC;
	t->size=1;
	t->treeNode[t->root].data=t->root;
	t->treeNode[t->root].dataCount_Sim=DATA_COUNT_EACH_NODE;
	for(int i=0;i<MAX_TREE_SIZE;i++){
		t->treeNode[i].father=NULL;
		t->treeNode[i].firstchild=NULL;
		t->treeNode[i].nextsibling=NULL;
		t->treeNode[i].dataCount_Sim=DATA_COUNT_EACH_NODE;
		t->treeNode[i].dataCount_Cal=DATA_COUNT_EACH_NODE;
		t->treeNode[i].energy=0;
		
	}
}
void checkTree(Tree *t){
	if(t->treeNode[2].firstchild->father==&t->treeNode[2]){
		cout<<"father-child right"<<endl;
	}else{
		cout<<"father-child wrong!!"<<endl;
	}
	if(t->treeNode[2].nextsibling->father==t->treeNode[2].father){
		cout<<"father-sibling right"<<endl;
	}else{
		cout<<"father-sibling wrong!!"<<endl;
	}
}
/********************End of creating a tree***********************************/
/********************E-RS Cal send**************************************/
double ps_ARQ(int sigma,double pk ){
	return 1-pow((1-pk),sigma+1);
	/*
	int i=0;
	double ps=0;
	for(i=1;i<=sigma+1;i++){
		ps+=pow((1-pk),i-1)*pk;	
	}
	return ps;
	*/
}
double ET_ARQ(int sigma,double pk,double bitSize,double dk){
	int i=0;
	double ps=0;
	double energy=0;
	ps=ps_ARQ(sigma,pk);
	/*
	if(pk==0.1){
		cout<<"sigma="<<sigma<<endl;
		cout<<"bitSize="<<bitSize<<endl;
		cout<<"dk"<<dk<<endl;
	}*/
//	for(i=1;i<=sigma+1;i++){
//			energy+=(i*pow((1-pk),i-1)*pk*(Etx(bitSize,dk)+Erx(bitSize)));
			/*
			if(pk==0.1){
				cout<<" pow((1-pk),i-1)="<<(i*pow((1-pk),i-1)*pk*(Etx(bitSize,dk)+Erx(bitSize)))<<endl;
			}*/
			
//	}
//	if(pk==0.1)cout<<endl;
	energy+=(1-(1+(sigma+1)*pk)*(pow((1-pk),(sigma+1))))*(Etx(bitSize,dk)+Erx(bitSize))/pk;
	energy+=ps*(Etx(ACK_SIZE*8,dk)+Erx(ACK_SIZE*8));
	energy+=(1-ps)*((sigma+1)*(Etx(bitSize,dk)+Erx(bitSize)));
	return energy;
}
double ps_RS(int sendNum,int M,double pk){
	int i=0;
	double ps=0;
	if(sendNum==0){
		cout<<"need debug *******************"<<endl;	
		return 1;
	}
	for(i=sendNum;i<=M;i++){
		ps+=CC(i-1,sendNum-1)*pow(pk,sendNum)*pow((1-pk),(i-sendNum));
	}
//debug
	if(ps>1.000001){
		cout<<"debug!!!!"<<"in ps_RS "<<"ps="<<ps<<" sendNum="<<sendNum<<" M=="<<M<<" pk=="<<pk<<endl;
	}
	return ps;
}
double AF_NC(double lastOneAlign,double pk,int sendNum, int M){
	int i,j;
	double af=0;
	double tmp_ps_RS=0;
	double tmp_for=0;
	double tmp_for_i=0;
	double tmp_for_j=0;
	af+=lastOneAlign*ps_ARQ(MAX_RETRY,pk);
	if(sendNum==0){
		
	}else if(sendNum==1){
		af+=MAX_APPEND_COUNT*ps_ARQ(MAX_RETRY,pk);
	}else{
		tmp_ps_RS=ps_RS(sendNum,M,pk);
		af+=MAX_APPEND_COUNT*sendNum*tmp_ps_RS;
		/*	
		for(i=1;i<=sendNum;i++){
			tmp_for+=i*CC(sendNum,i)*pow(pk,i)*pow((1-pk),(sendNum-i));
		}
		af+=MAX_APPEND_COUNT*(1-tmp_ps_RS)*tmp_for;	
		*/
		tmp_for_i=0;
		for(i=1;i<=sendNum-1;i++){
			tmp_for_j=0;
			for(j=0;j<=sendNum-i-1;j++){
				tmp_for_j+=CC((LAMBDA-1)*sendNum,j)*pow(pk,j)*pow((1-pk),((LAMBDA-1)*sendNum-j));
			}
			tmp_for_i+=i*CC(sendNum,i)*pow(pk,i)*pow((1-pk),(sendNum-i))*tmp_for_j;
		}
		af+=MAX_APPEND_COUNT*tmp_for_i;
	}
	return af;
}
double ET_RS(double lastOneAlign, int sendNum,int M,double pk,double bitSize,double dk){
	int i=0;
	double energy=0;
	double tmp_ps_RS=0;
	double tmp_for=0;
	double lastOneBitSize=0;
	

	if(sendNum==0){
		
	}else if(sendNum==1){
		energy+=ET_ARQ(MAX_RETRY,pk,bitSize,dk);
	}else{
		tmp_ps_RS=ps_RS(sendNum,M,pk);
		energy+=tmp_ps_RS*(Etx(ACK_SIZE*8,dk)+Erx(ACK_SIZE*8));
		tmp_for=0;
		for(i=sendNum;i<=M;i++){
			tmp_for+=(Etx(bitSize,dk)+Erx(bitSize))*i*CC(i-1,sendNum-1)*pow(pk,sendNum)*pow((1-pk),(i-sendNum));
		}
		energy+=tmp_for;
		tmp_for=0;
		energy+=(1-tmp_ps_RS)*M*(Etx(bitSize,dk)+Erx(bitSize));
	}
	//energy of last one
	if(lastOneAlign){
	//last one ARQ sending energy
		//formate a normal append frame
	
		lastOneBitSize=(25+lastOneAlign*SENSOR_DATA_SIZE)*8;
		energy+=ET_ARQ(MAX_RETRY,pk,lastOneBitSize,dk);
	}
	return energy;
}
void E_RS_send_cal(CSNode *p){
	//NC, append
	double arrivalDataCount=0;
	//after append
	int dataCountIntPart;
	//mind when sendNum=0
	int sendNum=0; //in sensor data block
	double lastOneAlign=0; //in float sensor data,it is a average float number
	//after NC
	int M=0;
	double energy=0;
	double bitSize=(25+1+1+1+MAX_APPEND_COUNT*SENSOR_DATA_SIZE)*8;
	double lastOneBitSize=(25+lastOneAlign*SENSOR_DATA_SIZE)*8;

	//check if is root
	if(p->father==NULL){
		return;
	}
	//append
	dataCountIntPart=p->dataCount_Cal/1;
	sendNum=(int)(dataCountIntPart/MAX_APPEND_COUNT);
	lastOneAlign=(dataCountIntPart)%MAX_APPEND_COUNT+(p->dataCount_Cal-dataCountIntPart);
	//end of append


	//coding
	M=LAMBDA*sendNum;
	//formate a NC frame
	bitSize=(25+1+1+1+MAX_APPEND_COUNT*SENSOR_DATA_SIZE)*8;

	//NC sending energy

	energy+=ET_RS(lastOneAlign,sendNum,M,PK_Cal,bitSize,p->dis_to_father);
	
	arrivalDataCount+=AF_NC(lastOneAlign,PK_Cal,sendNum,M);

	if(p->father==NULL){
	//this is root
	}else{
		p->father->energy+=(p->energy+energy);
		p->father->dataCount_Cal+=arrivalDataCount;
	}
}


/********************End of Cal send**************************************/
/********************E-RS Sim send**************************************/
void CTP_send_sim(CSNode *p){
	//ARQ,NOT append
	int arrivaldataCount_Sim=0;
	int randP;
	int reTryCount=MAX_RETRY+1;
	double energy=0;
	double bitSize=(25+SENSOR_DATA_SIZE)*8;
	//check if is root
	if(p->father==NULL)return ;
	for(int i=0;i<(p->dataCount_Sim);i++){
	//ARQ send each data
		while(reTryCount){
			reTryCount--;
			randP=random(PERCENT);
			energy+=Etx(bitSize,p->dis_to_father);
			energy+=Erx(bitSize);
			if(randP<PK_Sim){
			//delivered
				energy+=Etx(ACK_SIZE*8,p->dis_to_father);
				energy+=Erx(ACK_SIZE*8);
				arrivaldataCount_Sim++;
				reTryCount=0;
			}
		}
		reTryCount=MAX_RETRY+1;
	}
	if(p->father==NULL){
	//this is root
	}else{
		p->father->dataCount_Sim+=arrivaldataCount_Sim;
		arrivaldataCount_Sim=0;
		p->father->energy+=energy+p->energy;
		energy=0;
	}
}
void E_RS_send_sim(CSNode *p){

	int arrivaldataCount_Sim=0;
	int randP;
	int reTryCount=MAX_RETRY+1;
	bool sendDone=false;
	//fragment
	int sendNum=0;		//s
	int sendCount=0;	
	int lastOneAlign=0;	//send last data with ARQ
	int systemCode_rx=0;	//system code receive count

	//after E-RS
	int M=0;
	int tmp=0;
	int rx_count=0;
	double energy=0;
	double bitSize=(25+1+1+1+MAX_APPEND_COUNT*SENSOR_DATA_SIZE)*8;
	double lastOneBitSize=0;
	//p is root
	if(p->father==NULL) return ;
	//append
	lastOneAlign=(p->dataCount_Sim)%MAX_APPEND_COUNT;
	sendNum=(int)((p->dataCount_Sim-lastOneAlign)/MAX_APPEND_COUNT);
	if(sendNum==0){
	//append data less than one segment 
	//no send
	
	}else if(sendNum==1){
	//ARQ
	
		reTryCount=MAX_RETRY+1;
		while(reTryCount){
			reTryCount--;
			randP=random(PERCENT);
			energy+=Etx(bitSize,p->dis_to_father);
			energy+=Erx(bitSize);
			if(randP<PK_Sim){
			//delivered
				energy+=Etx(ACK_SIZE*8,p->dis_to_father);
				energy+=Erx(ACK_SIZE*8);
				reTryCount=0;
				arrivaldataCount_Sim+=sendNum*MAX_APPEND_COUNT;
			}
		}				
	}else{
	//E-RS	
		//encode
		M=LAMBDA*sendNum;
		//sendiing
		while(M&&!sendDone){
			randP=random(PERCENT);
			sendCount++;
			energy+=Etx(bitSize,p->dis_to_father);
			energy+=Erx(bitSize);
			if(randP<PK_Sim){
				//delivered
				rx_count++;
				if(sendCount<=sendNum){
					systemCode_rx++;
				}
				if(rx_count==sendNum){
					sendDone=true;
				}
			}else{
				//lost
			}
			M--;
		}
		if(sendDone==true){
		// ALL decode SUCCESS
			//ACK energy
			energy+=Etx(ACK_SIZE*8,p->dis_to_father);
			energy+=Erx(ACK_SIZE*8);
			arrivaldataCount_Sim+=sendNum*MAX_APPEND_COUNT;
		}else{
		//decode FAIl but system code SUCCESS
		//no ACK
			arrivaldataCount_Sim+=systemCode_rx*MAX_APPEND_COUNT;
		}
		
	}//end of else(sendNum)
	//send last align
	if(lastOneAlign){
		
		lastOneBitSize=(25+lastOneAlign*SENSOR_DATA_SIZE)*8;
		reTryCount=MAX_RETRY+1;
		while(reTryCount){
			reTryCount--;
			randP=random(PERCENT);
			energy+=Etx(lastOneBitSize,p->dis_to_father);
			energy+=Erx(lastOneBitSize);
			if(randP<PK_Sim){
			//delivered
				energy+=Etx(ACK_SIZE*8,p->dis_to_father);
				energy+=Erx(ACK_SIZE*8);
				arrivaldataCount_Sim+=lastOneAlign;
				reTryCount=0;
			}else{
			//lost
			}
		}
	}
	if(p->father==NULL){
	//p is root
	}else{
		p->father->dataCount_Sim+=arrivaldataCount_Sim;
		p->father->energy+=energy+p->energy;
	}
}
/********************End of Sim send**************************************/
void BFS(CSNode *n){
	CSNode *p;
	
	if(n==NULL)return ;
	if(n->firstchild!=NULL){
		BFS(n->firstchild);
		p=n->firstchild;
		while(p->nextsibling){
			BFS(p->nextsibling);
			p=p->nextsibling;
		}
	}
	
	//NOTE:should clear trace
	switch(choice){
		case SEND_CHOICE_CAL:{
			E_RS_send_cal(n);
			   }break;
		case SEND_CHOICE_SIM:{
			E_RS_send_sim(n);
			   }break;
		default :{}
	}
}
void clearTreeTrace(Tree *t){
	t->root=0;
	t->treeNode[t->root].dataCount_Sim=DATA_COUNT_EACH_NODE;
	t->treeNode[t->root].dataCount_Cal=DATA_COUNT_EACH_NODE;
	for(int i=0;i<MAX_TREE_SIZE;i++){
		t->treeNode[i].dataCount_Sim=DATA_COUNT_EACH_NODE;
		t->treeNode[i].dataCount_Cal=DATA_COUNT_EACH_NODE;
		t->treeNode[i].energy=0;
	}
}


int main(int argc, char **argv){
//0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95
//18
int i,j;
int pk_count=18;
int step;
double tmp_e=0;
double tmp_p=0;
double tmp_ar=0;
//result

double p_vs_pk[2][18]={0};
double e_vs_pk[2][18]={0};

double p_er[18]={0};
double e_er[18]={0};

	srand((unsigned)time(NULL));	

	for(i=0;i<RE_EXP;i++){
		initTree(&tree);
		creatTree(&tree,MAX_TREE_SIZE);
		for(step=0;step<pk_count;step++){
			//start of sim
			clearTreeTrace(&tree);
			choice=SEND_CHOICE_SIM;
			PK_Sim=(0.1+step*0.05)*PERCENT;
			tmp_e=0;
			tmp_p=0;
			tmp_ar=0;
			for(j=0;j<RE_RANDOM;j++){
					clearTreeTrace(&tree);
					BFS(&(tree.treeNode[tree.root]));	
					tmp_p+=(double)(tree.treeNode[tree.root]).dataCount_Sim/(double)MAX_TREE_SIZE;
					tmp_ar+=(tree.treeNode[tree.root]).dataCount_Sim;
					tmp_e+=(tree.treeNode[tree.root]).energy;
			}
			
			//average
			//start of sim
			tmp_e=tmp_e/(double)RE_RANDOM;
			tmp_ar=tmp_ar/(double)RE_RANDOM;
	//		tmp_e=tmp_e/tmp_ar;
			p_vs_pk[SEND_CHOICE_SIM-1][step]+=tmp_p/(double)RE_RANDOM;
			e_vs_pk[SEND_CHOICE_SIM-1][step]+=tmp_e;
			//end of sim
			
			//start of cal
			PK_Cal=(0.1+step*0.05);
			choice=SEND_CHOICE_CAL;
			clearTreeTrace(&tree);
			BFS(&(tree.treeNode[tree.root]));	
			p_vs_pk[SEND_CHOICE_CAL-1][step]+=(double)(tree.treeNode[tree.root]).dataCount_Cal/(double)MAX_TREE_SIZE;
//			e_vs_pk[SEND_CHOICE_CAL-1][step]+=(tree.treeNode[tree.root]).energy/(double)(tree.treeNode[tree.root]).dataCount_Cal;
			e_vs_pk[SEND_CHOICE_CAL-1][step]+=(tree.treeNode[tree.root]).energy;
			//end of cal
		}
	}
	//average Sim
	for(step=0;step<pk_count;step++){
		p_vs_pk[SEND_CHOICE_SIM-1][step]=p_vs_pk[SEND_CHOICE_SIM-1][step]/(double)(RE_EXP);
		e_vs_pk[SEND_CHOICE_SIM-1][step]=e_vs_pk[SEND_CHOICE_SIM-1][step]/(double)(RE_EXP);
		e_vs_pk[SEND_CHOICE_SIM-1][step]=e_vs_pk[SEND_CHOICE_SIM-1][step]/(double)1000; //nj to uj
	}
	//average Cal
	for(step=0;step<pk_count;step++){
		p_vs_pk[SEND_CHOICE_CAL-1][step]=p_vs_pk[SEND_CHOICE_CAL-1][step]/(double)(RE_EXP);
		e_vs_pk[SEND_CHOICE_CAL-1][step]=e_vs_pk[SEND_CHOICE_CAL-1][step]/(double)(RE_EXP);
		e_vs_pk[SEND_CHOICE_CAL-1][step]=e_vs_pk[SEND_CHOICE_CAL-1][step]/(double)1000; //nj to uj
	}
	//cal error
	for(step=0;step<pk_count;step++){
		//p
		tmp_p=0;
		if(p_vs_pk[SEND_CHOICE_CAL-1][step]>p_vs_pk[SEND_CHOICE_SIM-1][step]){
			tmp_p=p_vs_pk[SEND_CHOICE_CAL-1][step]-p_vs_pk[SEND_CHOICE_SIM-1][step];
		}else{
			tmp_p=p_vs_pk[SEND_CHOICE_SIM-1][step]-p_vs_pk[SEND_CHOICE_CAL-1][step];
		}
		p_er[step]=tmp_p/p_vs_pk[SEND_CHOICE_CAL-1][step];

		//e
		tmp_e=0;
		if(e_vs_pk[SEND_CHOICE_CAL-1][step]>e_vs_pk[SEND_CHOICE_SIM-1][step]){
			tmp_e=e_vs_pk[SEND_CHOICE_CAL-1][step]-e_vs_pk[SEND_CHOICE_SIM-1][step];
		}else{
			tmp_e=e_vs_pk[SEND_CHOICE_SIM-1][step]-e_vs_pk[SEND_CHOICE_CAL-1][step];
		}
		e_er[step]=tmp_e/e_vs_pk[SEND_CHOICE_CAL-1][step];
	}

	cout<<"simulation "<<RE_EXP*RE_RANDOM<<" Times"<<endl;
	cout<<"r=8"<<endl;
	cout<<"cita="<<SENSOR_DATA_SIZE<<endl;
	cout<<"x="<<MAX_APPEND_COUNT<<endl;
	cout<<"LAMBDA="<<LAMBDA<<endl;
	cout<<"Sigma="<<MAX_RETRY<<endl;
	//ouput
	cout<<"clear"<<endl;
//out put x
	cout<<"x=[";
	for(i=0;i<pk_count;i++){
		cout<<" "<<0.1+i*0.05;
	}
	cout<<"]"<<endl;
//out put p
	cout<<"p_RS_sim=[";
	for(i=0;i<pk_count;i++){
		cout<<setw(5)<<setprecision(2)<<setiosflags(ios::fixed)<<p_vs_pk[SEND_CHOICE_SIM-1][i];
	}
	cout<<"]"<<endl;
	cout<<"p_RS_cal=[";
	for(i=0;i<pk_count;i++){
		cout<<setw(5)<<setprecision(2)<<setiosflags(ios::fixed)<<p_vs_pk[SEND_CHOICE_CAL-1][i];
	}
	cout<<"]"<<endl;
//out put e
	cout<<"e_RS_sim=[";
	for(i=0;i<pk_count;i++){
		cout<<setw(10)<<setprecision(2)<<setiosflags(ios::fixed)<<e_vs_pk[SEND_CHOICE_SIM-1][i];
	}
	cout<<"]"<<endl;
	cout<<"e_RS_cal=[";
	for(i=0;i<pk_count;i++){
		cout<<setw(10)<<setprecision(2)<<setiosflags(ios::fixed)<<e_vs_pk[SEND_CHOICE_CAL-1][i];
	}
	cout<<"]"<<endl;
//out put error
	cout<<"p_err_RS=[";
	for(i=0;i<pk_count;i++){
		cout<<setw(8)<<setprecision(5)<<setiosflags(ios::fixed)<<p_er[i]*100;
	}
	cout<<"]"<<endl;
	cout<<"e_err_RS=[";
	for(i=0;i<pk_count;i++){
		cout<<setw(8)<<setprecision(5)<<setiosflags(ios::fixed)<<e_er[i]*100;
	}
	cout<<"]"<<endl;
	return 0;
}
