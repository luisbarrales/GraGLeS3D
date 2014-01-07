/* 
 * File:   main.cpp
 * Author: kuehbach
 *
 * Created on 27. März 2013, 07:52
 */

#include <cstdlib>
//#include <stdio.h>
#include <time.h>
//#include <string.h>
//#include <iostream>

#include "mymath.h"
#include "random.h"
#include "io.h"

#define MCRIT 15.0
#define TOL 0.1


using namespace std;

//generate n newOrientations from a reference according to a scatter
void randomSetOfOrientations (long n) {
	randomClass random;
	random.init ( -46356 );

	mathMethods mymath;

	cout << "oribucket_0.ori into CORe6NOO.ori 2 IMMOri.ori Converter at simulation time 0.0000" << endl;
	cout << "PHI2 22" << endl;
	cout << "9 0 0" << endl;
		
	double dvol = pow((double) (n), -1.0);
	double scatter = 3.489;
	for (long i = 0; i < n; i++) {
		double eulA[3];
		mymath.randomOrientation( eulA );
		//mymath.randomOri( eulA );
		//mymath.randomOriShoemake( eulA );
		
		cout << eulA[0] << "\t" << eulA[1] << "\t" << eulA[2] << "\t" << dvol << "\t" << scatter << endl;
	}
//	cout << "Last ID analyzed\n" << endl;
}


void eul2quat2eul_luis (long n) {
	//generate n random orientations, convert them into their quaternions and in Euler angle again
	//calculate occuring misorientation during each conversion step with misorientationCubic, ..quaternionCubic, 

	randomClass random;
	random.init ( -46356 );

	mathMethods mymath;

	cout << "Euler2Quaternion2EulerWithDisorientations_Luis testing" << endl;
	cout << "ID;ID;phi1;PHI;phi2;q0;q1;q2;q3;phi1;PHI;phi2;q0B;q1B;q2B;q3B;misorientationCubicEulerAB;misoriQxQAB;misoriQuatCubic0;misoriQuatCubic1;misoriQuatCubic2;misoriQuatCubic3" << endl;

	for (long i = 0; i < n; i++) {
		double eulA[3];
		double qA[4];
		double eulB[3];
		double qB[4];
		
		double misoriCubicEulAB;
		double quatDisori[4];
		double misoriQxQAB;

		mymath.randomOrientation( eulA );
		mymath.euler2quaternion( eulA, qA );
		mymath.quaternion2Euler( qA, eulB );
		mymath.euler2quaternion( eulB, qB );

		misoriCubicEulAB = mymath.misorientationCubic ( eulA[0], eulA[1], eulA[2], eulB[0], eulB[1], eulB[2] );
		mymath.misorientationQuaternionCubic ( qA, qB, quatDisori );
		misoriQxQAB = mymath.misorientationCubicQxQ( qA[0], qA[1], qA[2], qA[3],   qB[0], qB[1], qB[2], qB[3] );
		
		cout << "ID;" << i << ";" << eulA[0] << ";" << eulA[1] << ";" << eulA[2] << ";" << qA[0] << ";" << qA[1] << ";" << qA[2] << ";" << qA[3] << ";" << eulB[0] << ";" << eulB[1] << ";" << eulB[2] << ";"<< qB[0] << ";" << qB[1] << ";" << qB[2] << ";" << qB[3] << ";" << misoriCubicEulAB << ";" << misoriQxQAB << ";" << quatDisori[0] << ";" << quatDisori[1] << ";" << quatDisori[2] <<";" << quatDisori[3] << endl; 
	}
	cout << "Last ID analyzed\n" << endl;
}

void eulereuler2disori_luis (long n) {
	//generate n pairs of two random Euler angles, calculate disorientation between them
	randomClass random;
	random.init ( -46356 );

	mathMethods mymath;

	cout << "EulerEuler2Disori_Luis testing" << endl;
	//cout << "ID;ID;phi1A;PHIA;phi2A;phi1B;PHIB;phi2B;q0A;q1A;q2A;q3A;q0B;q1B;q2B;q3B;misorientationCubicEulerAB;misoriCubicEulerABOrigInv;misoriCubicEulerABCOReV2;misoriCubicCOREv2CorrectInvet;misoriQxQAB;misoriQuatCubic0;misoriQuatCubic1;misoriQuatCubic2;misoriQuatCubic3" << endl;
	cout << "ID;ID;phi1A;PHIA;phi2A;phi1B;PHIB;phi2B;q0A;q1A;q2A;q3A;q0B;q1B;q2B;q3B;misorientationCubicEulerAB;misoriQxQAB;misoriQuatCubic0;misoriQuatCubic1;misoriQuatCubic2;misoriQuatCubic3" << endl;

	for (long i = 0; i < n; i++) {
		double eulA[3];
		double eulB[3];
		double qA[4];
		double qB[4];
		
		double misoriCubicEulAB;
		//double misoriCubicEulABOrigInv;
		double misoriCubicCOReV2;
		//double misoriCubicCOReV2CorrectInvert;
		double quatDisori[4];
		double misoriQxQAB;

		mymath.randomOrientation( eulA );
		mymath.randomOrientation( eulB );
		
		mymath.euler2quaternion( eulA, qA );
		mymath.euler2quaternion( eulB, qB );

		
		misoriCubicEulAB = mymath.misorientationCubic ( eulA[0], eulA[1], eulA[2], eulB[0], eulB[1], eulB[2] );
		
		//misoriCubicEulABOrigInv = mymath.misorientationCubicOrigInv( eulA[0], eulA[1], eulA[2], eulB[0], eulB[1], eulB[2] );
		
		misoriCubicCOReV2 = mymath.misorientationCubicCOReV2( eulA[0], eulA[1], eulA[2], eulB[0], eulB[1], eulB[2] ); 
		
		//misoriCubicCOReV2CorrectInvert = mymath.misorientationCubicCorrectCOReV2InvertAndMult( eulA[0], eulA[1], eulA[2], eulB[0], eulB[1], eulB[2] ); 
		
		mymath.misorientationQuaternionCubic ( qA, qB, quatDisori );
	
		misoriQxQAB = mymath.misorientationCubicQxQ( qA[0], qA[1], qA[2], qA[3],   qB[0], qB[1], qB[2], qB[3] );
		
		//cout << "ID;" << i << ";" << eulA[0] << ";" << eulA[1] << ";" << eulA[2] << ";" << eulB[0] << ";" << eulB[1] << ";" << eulB[2] << ";" << qA[0] << ";" << qA[1] << ";" << qA[2] << ";" << qA[3] << ";" << qB[0] << ";" << qB[1] << ";" << qB[2] << ";" << qB[3] << ";" << misoriCubicEulAB << ";" << misoriCubicEulABOrigInv << ";" << misoriCubicCOReV2 << ";"  << misoriCubicCOReV2CorrectInvert << ";" << misoriQxQAB << ";" << quatDisori[0] << ";" << quatDisori[1] << ";" << quatDisori[2] <<";" << quatDisori[3] << endl;
		cout << "ID;" << i << ";" << eulA[0] << ";" << eulA[1] << ";" << eulA[2] << ";" << eulB[0] << ";" << eulB[1] << ";" << eulB[2] << ";" << qA[0] << ";" << qA[1] << ";" << qA[2] << ";" << qA[3] << ";" << qB[0] << ";" << qB[1] << ";" << qB[2] << ";" << qB[3] << ";" << misoriCubicEulAB << ";" << misoriQxQAB << ";" << quatDisori[0] << ";" << quatDisori[1] << ";" << quatDisori[2] <<";" << quatDisori[3] << endl;
	}
	cout << "Last ID analyzed\n" << endl;
}

void rotateBaby (long n) {
	//generates random Euler angle rotates that about random u,v,w axis about random angle
	randomClass random;
	random.init ( -46356 );

	mathMethods mymath;

	cout << "rotateBaby testing" << endl;
	cout << "ID;ID;phi1A;PHIA;phi2A;u;v;w;randomAngle;eulRot[0];eulRot[1];eulRot[2]" << endl;
	
	for (long i = 0; i < n; i++) {
		double eulA[3];
		double uvw[3];
		double angle;
		double eulRot[3];
		
		mymath.randomOrientation( eulA );
		double uvwnorm;
		do {
			uvw[0] = random.parkMiller();
			uvw[1] = random.parkMiller();
			uvw[2] = random.parkMiller();
			
			uvwnorm = sqrt(SQR(uvw[0]) + SQR(uvw[1]) + SQR(uvw[2]));
		} while ( uvwnorm <= 0.0 );

		//normalize
		uvw[0] /= uvwnorm;
		uvw[1] /= uvwnorm;
		uvw[2] /= uvwnorm;
			
		angle = random.parkMiller() * 2 * _PI_;
		do {
			angle = random.parkMiller() * 2 * _PI_;
		} while ( angle > 1.09606677); //angle of disorientation, which is in m-3m fundamental zone 62.8/360*2*_PI_

		mymath.rotateOrientation( eulA, angle, uvw[0], uvw[1], uvw[2], eulRot );
		
		cout << "ID;" << i << ";" << eulA[0] << ";" << eulA[1] << ";" << eulA[2] << ";" << uvw[0] << ";" << uvw[1] << ";" << uvw[2] << ";" << angle << ";" << eulRot[0] << ";" << eulRot[1] << ";" << eulRot[2] << endl; 
	}
	cout << "Last ID analyzed\n" << endl;
}

void rotateIntoConeBaby (long n) {
	//generates random Euler angle rotates that about random u,v,w axis about random angle <= max angle
	randomClass random;
	random.init ( -46356 );

	mathMethods mymath;

	cout << "rotateBabyIntoFixedCone testing" << endl;
	cout << "ID;ID;phi1A;PHIA;phi2A;u;v;w;randomAngle;eulRot[0];eulRot[1];eulRot[2];maxAngle;disoriCubicQxQ" << endl;
	
	for (long i = 0; i < n; i++) {
		double eulA[3];
		double uvw[3];
		double angle;
		double eulRot[3];
		
		mymath.randomOrientation( eulA );
		uvw[0] = 1;
		uvw[1] = 1;
		uvw[2] = 2;
		angle = random.parkMiller() * 2 * _PI_;

		double maxAngle =  0.261799387; //15degrees
		cout << maxAngle << endl;

		mymath.newOrientationFromReferenceFixedAngularCone( eulA, maxAngle , angle, uvw[0], uvw[1], uvw[2], eulRot );
		
		double qA[4];
		double qB[4];
		double disoriQxQ;
		mymath.euler2quaternion( eulA, qA );
		mymath.euler2quaternion( eulRot, qB );
		disoriQxQ = mymath.misorientationCubicQxQ( qA[0], qA[1], qA[2], qA[3], qB[0], qB[1], qB[2], qB[3] );
		
		cout << "ID;" << i << ";" << eulA[0] << ";" << eulA[1] << ";" << eulA[2] << ";" << uvw[0] << ";" << uvw[1] << ";" << uvw[2] << ";" << angle << ";" << eulRot[0] << ";" << eulRot[1] << ";" << eulRot[2] << ";" << maxAngle << ";" << disoriQxQ << endl; 
	}
	cout << "Last ID analyzed\n" << endl;

}

void iamAquatDistance (long n)
{
	//generates two random euler angles, converts into quaternions and calculates disori
	randomClass random;
	random.init ( -46356 );

	mathMethods mymath;

	cout << "iamAquatDistance testing" << endl;
	cout << "ID;ID;phi1A;PHIA;phi2A;phi1B;PHIB;phi2B;q0A;q1A;q2A;q3A;q0B;q1B;q2B;q3B;distanceAB" << endl;
	
	for (long i = 0; i < n; i++) {
		double eulA[3];
		double eulB[3];
		double qA[4];
		double qB[4];
		double distanceAB;

		mymath.randomOrientation( eulA );
		mymath.randomOrientation( eulB );
		mymath.euler2quaternion( eulA, qA );
		mymath.euler2quaternion( eulB, qB );

		distanceAB = mymath.distanceBetweenQuaternions( qA, qB );
		
		cout << "ID;" << i << ";" << eulA[0] << ";" << eulA[1] << ";" << eulA[2] << ";" << eulB[0] << ";" << eulB[1] << ";" << eulB[2] << ";" << qA[0] << ";" << qA[1] << ";" << qA[2] << ";" << qA[3] << ";" << qB[0] << ";" << qB[1] << ";" << qB[2] << ";" << qB[3] << ";" << distanceAB << endl; 
	}
	cout << "Last ID analyzed\n" << endl;
}

void quat2rodrigues (long n)
{
	//generates random euler angles, converts into quaternion and calculates Rodrigues type angle axis representation
	randomClass random;
	random.init ( -46356 );

	mathMethods mymath;

	cout << "quat2rodrigues testing" << endl;
	cout << "ID;ID;phi1A;PHIA;phi2A;q0A;q1A;q2A;q3A;angle;qr1;qr2;qr3" << endl;
	
	for (long i = 0; i < n; i++) {
		double eulA[3];
		double qA[4];
		double qr[4];

		mymath.randomOrientation( eulA );
		mymath.euler2quaternion( eulA, qA );
		mymath.quaternion2rodrigues( qA, qr );
		
		cout << "ID;" << i << ";" << eulA[0] << ";" << eulA[1] << ";" << eulA[2] << ";" << qA[0] << ";" << qA[1] << ";" << qA[2] << ";" << qA[3] << ";" << qr[0] << ";" << qr[1] << ";" << qr[2] << ";" << qr[3] << endl; 
	}
	cout << "Last ID analyzed\n" << endl;
}


void eulscatter (long n)
{
	//generates random euler angle and random scatter around that to form new orientation
	randomClass random;
	random.init ( -46356 );

	mathMethods mymath;

	cout << "eulerscatter testing" << endl;
	cout << "ID;ID;phi1A;PHIA;phi2A;minDev;myDev;eulNew0;eulNew1;eulNew2;disoriANew" << endl;
	
	for (long i = 0; i < n; i++) {
		double eulA[3];
		//mymath.randomOrientation( eulA );
		mymath.randomOriShoemake( eulA );
		double eulNew[3];
		
		double min_myDev = 0.017453292; //1degree is a resonable lower limit
		double myDev = random.parkMiller() * 2 * _PI_;
		
		while (myDev < min_myDev) {
			myDev = random.parkMiller() * 2 * _PI_;
		}

		mymath.newOrientationFromReference(eulA, myDev, eulNew );
		
		double disoriANew = mymath.misorientationCubic( eulA[0], eulA[1], eulA[2], eulNew[0], eulNew[1], eulNew[2] );
		
		cout << "ID;" << i << ";" << eulA[0] << ";" << eulA[1] << ";" << eulA[2] << ";" << min_myDev << ";" << myDev << ";" << eulNew[0] << ";" << eulNew[1] << ";" << eulNew[2] << ";" << disoriANew << endl; 
	}
	cout << "Last ID analyzed\n" << endl;
}

void bubbletest ( void )
{
	mathMethods mymath;

	double test[4] = {3.0, 4.0, 1.0, 2.0};

	mymath.bubbleSort(test , 4);

	cout << test[0] << ";" << test[1] << ";" << test[2] << ";" << test[3] << endl;
}


int main(int argc, char** argv) {
/*	randomClass random;
	random.init(-46356);
   
    double cube[3], plage[3], nd45[3], etest[3];
    double qc[4], qp[4], qnd[4], qtest[4], qrr[4];
    //Cube
    cube[0] = 0.0 / 180 * _PI_;
    cube[1] = 0.0 / 180 * _PI_;
    cube[2] = 0.0 / 180 * _PI_;
    //P
    plage[0] = 65.0 / 180 * _PI_;
    plage[1] = 45.0 / 180 * _PI_;
    plage[2] = 0.0 / 180 * _PI_;
    //45ND
    nd45[0] = 45.0 / 180 * _PI_;
    nd45[1] = 0.0 / 180 * _PI_;
    nd45[2] = 0.0 / 180 * _PI_;
    
    double m_cube, m_plage, m_nd45;
       
    ioP myio;
	
    dataBlockP ori = myio->readDataBlock("OriIMM", "oribucket.10.IMM.uds");
    
    int n = ori->lineCount;
    dataLineP line = ori->first;
    euler2quaternion(cube, qc);
    euler2quaternion(plage, qp);
    euler2quaternion(nd45, qnd);
    
    double sum_cube, sum_plage, sum_nd45;
    sum_cube = 0;
    sum_plage = 0;
    sum_nd45 = 0;
    
   cout << "Cube\t\tPLage\t\t45ND" << endl;    
   for (int i = 0; i < n - 1; i++) {   
       
        etest[0] = myio->getReal(line,1) / 180 * _PI_;
        etest[1] = myio->getReal(line,2) / 180 * _PI_;
        etest[2] = myio->getReal(line,3) / 180 * _PI_;
        euler2quaternion(etest, qtest);
        
        //check cube
        misorientationQuaternionCubic(qtest,qc, qrr);
        m_cube = acos(qrr[0])*2  * 180 / _PI_;
        if (m_cube < MCRIT)
            sum_cube += myio->getReal(line,4);
        
        
        //check P
        misorientationQuaternionCubic(qtest,qp, qrr);
        m_plage = acos(qrr[0])*2  * 180 / _PI_;
        if (m_plage < MCRIT)
            sum_plage += myio->getReal(line,4);
        
        //check 45ND
        misorientationQuaternionCubic(qtest,qnd, qrr);
        m_nd45 = acos(qrr[0])*2  * 180 / _PI_;
        if (m_nd45 < MCRIT)
            sum_nd45 += myio->getReal(line,4);
        
            
        
       //cout << m_cube << "\t\t;" << m_plage << "\t\t;" << m_nd45 << endl;
         
        // if (i % 20 == 0)
        //     cout << endl;
      
        line = line->next;    
    }
       
    cout << sum_cube << "\t\t" << sum_plage << "\t\t" << sum_nd45 << endl;
    
    //check point to point misorientation successively
    /*
    ioP myio;
    dataBlockP ori = myio->readDataBlock("OriIMM", "oribucket.10.IMM.uds");
    
    int n = ori->lineCount;
    dataLineP line = ori->first;
    
    double e1[3], e2[3];
    e1[0] = myio->getReal(line,1) / 180 * _PI_;
    e1[1] = myio->getReal(line,2) / 180 * _PI_;
    e1[2] = myio->getReal(line,3) / 180 * _PI_;
   // cout << e1[0] << "," << e1[1] << "," << e1[2] << endl;
    
   double q[4], p[4], qrr[4], misori;  //Euler angles are in rad
  
    
    for (int i = 0; i < n - 1; i++) {
        line = line->next;
        
        euler2quaternion(e1, q);
             
        e2[0] = myio->getReal(line,1) / 180 * _PI_;
        e2[1] = myio->getReal(line,2) / 180 * _PI_;
        e2[2] = myio->getReal(line,3) / 180 * _PI_;
        
        euler2quaternion(e2, p);
        
        misorientationQuaternionCubic(q,p, qrr);
       
         misori = acos(qrr[0])*2  * 180 / _PI_;

         cout << misori << ";"; 
         
         if (i % 20 == 0)
             cout << endl;
        
         e1[0] = myio->getReal(line,1) / 180 * _PI_;
        e1[1] = myio->getReal(line,2) / 180 * _PI_;
        e1[2] = myio->getReal(line,3) / 180 * _PI_;
     //   cout << e1[0] << "," << e1[1] << "," << e1[2] << endl;
        //line = line->next;    
    }
    
     * 
     * /
    //TEST CORE QUATERNION CALCULATIONS
    /*
    double eref[3];
    eref[0] = 290.0 / 180 * _PI_;
    eref[1] = 290.0 / 180 * _PI_;
    eref[2] = 290.0 / 180 * _PI_;
   
    cout << "triplett;phi1;PHI;phi2;eres1;eres2;eres3" << endl;
   
   double angle = 1.0;
   double misori = 0;
   
   int n = 100;
   for (int triplett = 1; triplett <= n; triplett++ ) {
        
        double eresult[3];
        
    
        newOrientationFromReference(eref, angle, eresult, &random);
       
        //double q[4], p[4], qrr[4];  //Euler angles are in rad
        //euler2quaternion(eref, q);
        //euler2quaternion(eresult, p);
       
        
        //misorientationQuaternionCubic(q,p, qrr);
        //cout << q[0] << ";" << q[1] << ";" << q[2] << ";" << q[3] << endl;
        //cout << p[0] << ";" << p[1] << ";" << p[2] << ";" << p[3] << endl;
        //double misori = 1; //acos(qrr[0])*2;
        double erefrad[3], q[4], p[4], qmis[4];
        erefrad[0] = eref[0];// / 180 * _PI_;
        erefrad[1] = eref[1];// / 180 * _PI_;
        erefrad[2] = eref[2];// / 180 * _PI_;
        
        euler2quaternion(erefrad, q);
        euler2quaternion(eresult,p);
        misorientationQuaternionCubic(q,p,qmis);
        misori = qmis[0];
        misori = acos(misori);
        misori *= 2 / _PI_ * 180;
        
                
       // misori = 0.0;
       // misori = misorientationCubic( eref[0], eref[1], eref[2], eresult[0], eresult[1], eresult[2] );
       //misori *= 180;
        //misori /= _PI_;
           
    cout << triplett << ";" << eref[0] << ";" << eref[1] << ";" << eref[2] << ";" << eresult[0] << ";" << eresult[1] << ";" << eresult[2] << ";"  << misori << endl;
            //<< qresult[0] << ";" << qresult[1] << ";" << qresult[2] << ";" << qresult[3] << ";"
   }
   
   
    
    //test randomMisorienation
    /*double qres[4];
    double angle = 4.0;
    double q0 = 0;
    int n = 100;
    for (int test = 1; test<=n; test++) {
        randomMisorientation(angle, qres, &random);
        q0 = acos(qres[0])*2*180/_PI_;
        cout << q0 << "\t" << qres[1] << "\t" << qres[2] << "\t" << qres[3] << endl; //<< "\t";
        //if (test % 10 == 0) cout << endl;
    }*/
	
	long n = 1;
	long mode = 1;
	
	n = atoi(argv[1]);
	mode = atoi(argv[2]);

	//bubbletest();
	cout << acos(0.0) << endl;
	cout << acos(-1.0) << endl;
	cout << acos(1.0) << endl;
	
	
	if (mode == 0) randomSetOfOrientations( n );
	
	if (mode == 1) eul2quat2eul_luis( n );

	if (mode == 2) eulereuler2disori_luis( n );

	if (mode == 3) rotateBaby( n );

	if (mode == 4) iamAquatDistance ( n );

	if (mode == 5) eulscatter ( n );

	//if (mode == 6) rotateIntoConeBaby( n ); //###MK::currently makes problems
	

	return 0;
}




//convert n Euler angles in quaternion and back in Euler angles
/*
int main(int argc, char** argv) {
    randomClass random;
    random.init(-46356);
    //138.89 , 0 , 0 , the quaternion newly generated0.997322 , 0.0659125 , 0.026475 , 0.0174099

   double e[3];
   double res[18];
   int n = 1000;
   //create n Euler triplets, convert them in a quaternion and convert them back
   cout << "triplett;phi1;PHI;phi2;q0;q1;q2;q3;phi1';PHI';phi2';q0';q1';q2';q3';qmis1;qmis2;qmis3;qmis4" << endl;
   
   for (int triplett = 1; triplett <= n; triplett++ ) {
       e[0] = 2 * _PI_ * random.parkMiller();
       e[1] = 1 * _PI_ * random.parkMiller();
       e[2] = 2 * _PI_ * random.parkMiller();
       
       double q[4];      
       eulerToQuaternion(e[0], e[1], e[2], q);
       
       res[0] = e[0];
       res[1] = e[1];
       res[2] = e[2];
        res[3] = q[0]; //log quaternion components temporarily
        res[4] = q[1];
        res[5] = q[2];
        res[6] = q[3];
        
        double ress[3];
        
        quaternion2Euler(q, ress);
        res[7] = ress[0];
        res[8] = ress[1];
        res[9] = ress[2];
        
        double p[3];
        double quatMisori[4];
        eulerToQuaternion(ress[0], ress[1], ress[2], p);
        res[10] = p[0]; //log quaternion components temporarily
        res[11] = p[1];
        res[12] = p[2];
        res[13] = p[3];
               
        misorientationQuaternionCubic(q,p,quatMisori);
        res[14] = quatMisori[0]; //log quaternion components temporarily
        res[15] = quatMisori[1];
        res[16] = quatMisori[2];
        res[17] = quatMisori[3];        
        
        
        //output
       cout << triplett << ";";
        for (int i = 0; i<18; i++)
            cout << res[i] << ";";
       cout << endl;
      
      
        //check if eul2quat2eul yields the same results within a tolerance
        /*if ( ( fabs(e[0] - ress[0]) < TOL ) && ( fabs(e[1] - ress[1]) < TOL ) && ( fabs(e[2] - ress[2]) < TOL ) ) 
            if ( ( (e[0] + ress[0]) - (e[2] + ress[2]) ) < TOL )
                cout << ";success" << endl;
            else
                cout << ";FAIL!" << endl;
        else
            cout << ";FAIL!" << endl;*/
//   }
/* 
   //new orientation from reference bogenmass
   
    //e[0] = 355.0;// / 180.0 * _PI_;
    //e[1] = 0.0;// / 180.0 * _PI_;
    //e[2] = 0.0;// / 180.0 * _PI_;
     
   
   
   double q[4];
   
   //for eulerToQuaternion
    e[0] = -12.0 / 180.0 * _PI_;
    e[1] = 180.0 / 180.0 * _PI_;
    e[2] = -12.0 / 180.0 * _PI_;
    
    
    
    double angle = 10;
    
    double q[7] = {-9999, -9999 , -9999, -9999, -9999 , -9999, -9999};
    
    //newOrientationFromReference( e, angle, q );
    
    eulerToQuaternion(e[0], e[1], e[2], q);
    
    double result[7] = {-9999, -9999 , -9999, -9999, -9999 , -9999, -9999};
    
        result[3] = q[0]; //log quaternion components temporarily
        result[4] = q[1];
        result[5] = q[2];
        result[6] = q[3];
        
    //qBungeZXZ = - qMatthies (MTEX)
    quaternion2Euler(q, result); //>(rad)
    
    cout << result[0] << " , " << result[1] << " , " << result[2] << endl;
    
    result[0] *= 180 / _PI_; // (rad) > (deg)
    result[1] *= 180 / _PI_;
    result[2] *= 180 / _PI_;
    
    
    
    //result[3] = q[3]; //transfer q to result
    //result[4] = q[4];
    //result[5] = q[5];
    //result[6] = q[6];
     
    
    
    
    cout << result[0] << " , " << result[1] << " , " << result[2] << " , the quaternion newly generated " << result[3] << " , " << result[4] << " , " << result[5] << " , " << result[6] << endl;

    
    return 0;
}
*/

