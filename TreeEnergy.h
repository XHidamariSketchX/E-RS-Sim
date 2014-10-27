/*
 * TreeEnergy.h
 *
 *  Created on: 2014Äê3ÔÂ2ÈÕ
 *      Author: xuji
 */

#ifndef TREEENERGY_H_
#define TREEENERGY_H_

#define MAX_TREE_SIZE 200
#define DEGREE 10

#define LOC_SIZE 1000 //1km*1km
#define ROOT_LOC LOC_SIZE/2 //root position (500,500)
#define DIS 100 //transmission distance
//#define MIN_DIS 30
#define DATA_COUNT_EACH_NODE 1
#define PERCENT 100
#define MAX_RETRY 5 	//sigma
#define LAMBDA 5			//lambda
#define MAX_APPEND_COUNT 10	//m segment size
#define SENSOR_DATA_SIZE 8 //cita in Bytes
//NOTE SENSOR_DARA_SIZE*MAX_APPEND_COUNT<=100

#define ACK_SIZE 5//5 Bytes

#define RE_EXP 1000
#define RE_RANDOM 100

#define SEND_CHOICE_CAL 1
#define SEND_CHOICE_SIM 2

typedef struct Node {
	int data;	//node id
	int dataCount_Sim;	// x SIM
	double dataCount_Cal;//x CAL
	double energy;
	int x;
	int y;
	double dis_to_father;	//dk
	struct Node * firstchild, *nextsibling,* father;

}CSNode;
typedef struct t{
	int size;
	int root;
	CSNode treeNode[MAX_TREE_SIZE];
}Tree;

#endif /* TREEENERGY_H_ */
