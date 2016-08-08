#include<iostream>
#include<fstream>
#include<math.h>

using namespace std;

const int rows1 = 152;
const int rows2 = 197;
const int rows_common = (rows1 + rows2);
const int cols = 2;

void common(double [rows1][cols],double [rows2][cols],double [rows_common][cols+1],int&);

void min(double [rows_common][cols+1],double [rows_common][cols+1],int);

void chang_miller(double [rows_common][cols+1],double [rows_common][cols], double, double, int ,int);

void oniom_xs(double [rows_common][cols+1],double [rows_common][cols],int ,int );

int main(int argc, char* argv[]){
  if(argc<2){
      cerr << "This is a function which can be used to get common points between three files" << endl;
  }
  else{
    ifstream myReadFile;
    myReadFile.open(argv[1]);
    int number_of_lines1 = 0;
      string line;
      ifstream myfile(argv[1]);
    while (getline(myfile, line))
      ++number_of_lines1;
    double output1, output2;
    double  points1 [number_of_lines1][2];
    if (myReadFile.is_open()) {
      for(int i=0; i < number_of_lines1; i++){
        myReadFile >> output1 >> output2;
        points1[i][0] = output1;
        points1[i][1] = output2;
      }
    }
    myReadFile.close();
    myReadFile.open(argv[2]); // reading in excluded topologies
    int number_of_lines2 = 0;
    string line2;
    ifstream myfile2(argv[2]);
    while (getline(myfile2, line2))
      ++number_of_lines2;
    double  points2 [number_of_lines2][2];

    if (myReadFile.is_open()) {
      for(int i=0; i < number_of_lines2; i++){
        myReadFile >> output1 >> output2;
        points2[i][0] = output1;
        points2[i][1] = output2;
      }
    }
    myReadFile.close();
 
    double common_array[rows_common][cols+1];
    int num = 0;
    common(points1,points2,common_array,num);

    double shifted_to_zero[rows_common][cols+1];
    min(common_array,shifted_to_zero,num);

//chang_miller   

    double chang_miller_energy[rows_common][cols];
    double A = 0.1;
    double alpha = 0.1;
    int autosub = 100;
    
    chang_miller(common_array,chang_miller_energy,A,alpha,autosub,num);

 //   Now ONIOM-XS
   double oniom_xs_energy[rows_common][cols];
   oniom_xs(common_array, oniom_xs_energy,autosub,num);
   for(int i=0;i<num;i++){
     cout << oniom_xs_energy[i][0] << " " << oniom_xs_energy[i][1]  << endl;
   }

  }
  return 0;
}

void oniom_xs(double common[rows_common][cols+1],double oniom_xs_energy[rows_common][cols],int autosub,int num){
  double x[num];
  double w[num];
    for(int i=0; i<num; i++){
      x[i] = (common[i][0] - common[autosub][1])/(common[i][1] - common[autosub][1]);
 
      w[i] = 6*(pow((x[i] - 0.5),5)) - 5*(pow((x[i] - 0.5),3)) + (15/8)*(x[i] - 0.5) + 0.5;
 
      oniom_xs_energy[i][0] = (1 - w[i])*common[i][0] + w[i]*common[i][1];
      oniom_xs_energy[i][1] = common[i][2];
    }
}

void chang_miller(double common[rows_common][cols+1],double chang_miller_energy[rows_common][cols], double A, double alpha,int autosub,int num){
  double square_non_diag_CM[rows_common];
  for(int i=0; i<num; i++){
    square_non_diag_CM[i] = A*exp(-1*alpha*(common[i][2]-common[autosub][2])*(common[i][2]-common[autosub][2]));
    chang_miller_energy[i][0] = 0.5*(common[i][0] + common[i][1] - sqrt((common[i][0]-common[i][1])*(common[i][0]-common[i][1]) + 4*square_non_diag_CM[i]));
    chang_miller_energy[i][1] = common[i][2];
  }
}

void min(double common[rows_common][cols+1],double shifted_to_zero[rows_common][cols+1],int number_of_common_points){
  double min[2];
  min[0] = common[0][0];
  min[1] = common[0][1];
  for(int i=0; i < number_of_common_points; i++){
    if(common[i][0] < min[0])
      min[0] = common[i][0];
    if(common[i][1] < min[1])
      min[1] = common[i][1];
  }
  for(int i=0; i < number_of_common_points; i++){
    if((common[i][0] != 0)&&(common[i][1] != 0)&&(common[i][2] != 0)){
      shifted_to_zero[i][0] = common[i][0] - min[0];
      shifted_to_zero[i][1] = common[i][1] - min[1];
    }
    shifted_to_zero[i][2] = common[i][2];
  }
}

void common(double points1[rows1][cols],double points2[rows2][cols],double common[rows_common][cols+1],int& number_of_common_points){
    double points[rows_common][cols];
    for(int i=0 ;i < rows1; i++){
        points[i][0] = points1[i][0];
        points[i][1] = points1[i][1];
    }

    for(int i=0; i< rows2; i++){
        points[i + rows2][0] = points2[i][0];
        points[i + rows2][1] = points2[i][1];
    }

  for (int i = 0; i < rows_common; i++){
    for (int j = i + 1;j < rows_common; j++){
      if ((points[i][1] == points[j][1])&&(points[i][1] != 0)){
        common[i][0] = points[i][0];
        common[i][1] = points[j][0];
        common[i][2] = points[i][1];
        number_of_common_points++;
      }
    }
  }
}
