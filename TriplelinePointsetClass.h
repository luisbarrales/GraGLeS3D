
/* 
 * File:   TriplelinePointsetClass.h
 * Author: Arbeit
 *
 * Created on 18. Dezember 2016
 */

#ifndef TRIPLELINEPOINTSETCLASS_H
#define TRIPLELINEPOINTSETCLASS_H

#include <vector>
#include <Eigen/Dense>
#include "Spline.h"

class Point3D;
class LocalPointset;

class InterfacialElement;

using namespace std;


class  TriplelinePointsetClass{
public:
    
    
    TriplelinePointsetClass();
    TriplelinePointsetClass(const Eigen::MatrixXd*, int);
    TriplelinePointsetClass(vector<Eigen::Vector3d>&, vector<InterfacialElement*>);
    ~TriplelinePointsetClass();
  
    void use_PointCloudGenerator(int);
    void input_TPS_OutOfTxtFile(string);
    void input_TPS_LevelSet(vector<Eigen::Vector3d>&, vector<InterfacialElement*>);
    void set_Parameter();
    void process_TriplelinePointset(); // const oder nicht const... leiber neue matrizen erstellen ??
    
    
    inline Eigen::MatrixXd* get_PointToPointDistancesMatrix(){
        return PointToPointDistances;
    }
   
    inline int get_N_TriplelinePoints(){
        return N;
    }

    inline vector<Point3D>* get_TriplelinePointsetBeforeMLS(){
        return TPS_BeforeMLS;
    }
   
    inline vector<Point3D>* get_TriplelinePointsetAfterMLS(){
        return TPS_AfterMLS;
    }
    
    
    inline int get_n(){
        return n;
    }
    
  
    inline int get_count(){
        return count;
    }
    
    inline int get_iterationcount(){
        return iterationscount;
    }
    inline double get_length(){
    	return lenght_spline;
    }
    inline vector<Point3D>* get_TPS_processed(){
    	return TPS_processed;
    }
    
    
private:
    
    //output functions:
    void output_internFunctionIntoTxtFile();
    void output_TPS_IntoTxtFile_onlyCoordinates(vector<Point3D>*, string, int);
    void output_TPS_IntoTxtFile(vector<Point3D>*, string, int);
    void output_EigenMatrixXd(Eigen::MatrixXd*, string);
    void output_LP_Object_IntoTxtFile(double, double);
    void output_DataDebugging_txt(string);
    
    //process_TriplelinePointset():
    //
    //MLS:
    void calc_PointToPointDistances(vector<Point3D>*);
    void smoothe_TPS();
    void calc_MLS_Iteration();
    void calc_Sufficient_LP_forMLS(int);
    void calc_NextNeighborSet_RadiusH(double, vector<Point3D>*);
    void switch_Buffer(); //nach einer iteration über das komplette Punktset. output des vorherigen schrittes als input des nächsten setzen
    
    //Ordering:
    void orderAndReduce_TPS();
        void calc_starting_point_ordering();
        int calc_RandomID(int);
        int calc_Next_Neighbor_ID(int);
        double calc_Next_Neighbor_Distance(int);
        int check_existing_neighbor(int, int);
        bool check_direction_dotProduct(int, int, int);
        int check_Neighbor_within_min_gap();
        void ordering_dynamicDistances_method();
        	void calc_Sufficient_LP_forOrdering(int);
            void order_direction_forward(int);
            void order_direction_backward(int);
            void copy_Buffer_ReducedAndOrdered();
        void ordering_Next_Neighbor_method();
            void calc_nonordered_Neighbors_forward(int);
            void calc_nonordered_Neighbors_backward(int);
    
    //Spline:
    void calc_BSpline_ProcessedTPS();  
    
    
    //variables+containers:
    //MovingLeastSquares:= MLS
    vector<Point3D>* TPS_BeforeMLS;
    vector<Point3D>* TPS_AfterMLS;
    vector<Point3D>* TPS_Ordered_Forward_Buffer;
    vector<Point3D>* TPS_Ordered_Backward_Buffer;
    vector<Point3D>* TPS_processed;
    
    
    Eigen::MatrixXd* BoolListPointsMovedInMLSIteration; //Information dass Punkt gemoved wurde speichern!!!! //vermerken ob punkt überhaupt berechnet wurde
    Eigen::MatrixXd* AverageApproxmationError;
    Eigen::MatrixXd* PointToPointDistances; //distances r=|deltaP|
    
    LocalPointset* LP_Object; //statisch alloc mit vorgegebener Vectorgröße
    splineclass::Spline* Spline_Object;
    
    // pointcloud generator 
    int PC_GEN; 
    
    //global variables for txt file output name declarations:
    ofstream* pOutFile;
    int count; //set in calc_SufficientLocalPointsetForLoopPoint
    int iterationscount; //set in calc_SufficientLocalPointsetForLoopPoint
  
    int n0; //vorgesehene Anzahl Schleifendurchläufe über das gesamte Triplelinepunktset
    int n; //aktueller Schleifendurchlauf
    double H0;
    int N; //Punktanzahl im TriplelinePointsetInput
    int maxPointID;
    int Nmin; //mindestanzahl Punkte in lokaler regression bzw local point set A
    double dH; // increment H
    double rho0;
    double epsilon0; //prescribed  average approximation error to quatratic regression line during moving least squares 
    
    //global variables ordering:
    double rho_order;
    double weight_order;
    double min_gap; 
    double max_gap;
    int ID_starting_point;
    int ID_last_point;
    int ID_ordered_forward;
    int ID_ordered_backward;
    
    //global variables spline:
    double lenght_generatorFct;
    double lenght_spline;
   
    //Rechenmatrizenm
    Eigen::Vector3d* Vector3d_1;
    Eigen::Vector3d* Vector3d_2;
    
 
};

#endif /* NEWCLASS_H */

