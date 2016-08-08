#include<iostream>
#include<fstream>
#include<math.h>

using namespace std;

const int rows1 = 183;
const int rows2 = 152;
const int rows3 = 183;
const int rows_common = (rows1 + rows2 + rows3);
const int cols = 2;

void common(double [rows1][cols],double [rows2][cols],double [rows3][cols],double [rows_common][cols+2],int&);

void min(double [rows_common][cols+2],double [rows_common][cols+2],int);

void shifting_to_autosub(double [rows_common][cols+2],double [rows_common][cols+2],double,int);

void evb(double [rows_common][cols+2],double [rows_common][cols+2],int);

int main(int argc, char* argv[]){
  if(argc<2){
      cerr << "This is a function which can be used to get evb energy. It takes three energy-distance files as input arguments" << endl;
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
    myfile.close();
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
    myfile2.close();
    myReadFile.close();
    myReadFile.open(argv[3]); // reading in trimer
    int number_of_lines3 = 0;
    string line3;
    ifstream myfile3(argv[3]);
    while (getline(myfile3, line3))
      ++number_of_lines3;
    double  points3 [number_of_lines3][2];
    if (myReadFile.is_open()) {
      for(int i=0; i < number_of_lines3; i++){
        myReadFile >> output1 >> output2;
        points3[i][0] = output1;
        points3[i][1] = output2;
      }
    }
    myReadFile.close();
    double common_array[rows_common][cols+2];
    int num = 0;
    common(points1,points2,points3,common_array,num);

    double shifted_to_zero[rows_common][cols+2];

    min(common_array,shifted_to_zero,num);

    ofstream myfile4;
    myfile4.unsetf (ios :: floatfield);
    myfile4.precision(15);
    myfile4.open("shifted_to_zero.txt"); 

    for(int i=0;i<num;i++){
      myfile4 << shifted_to_zero[i][0] << " " << shifted_to_zero[i][1] << " " << shifted_to_zero[i][2] << " " << shifted_to_zero[i][3] << endl;
    } 

    myfile4.close();

    double shifted_to_autosub[rows_common][cols+2];
    double autosub_distance = 1.386929;

    shifting_to_autosub(shifted_to_zero,shifted_to_autosub,autosub_distance,num);

    ofstream myfile5;
    myfile5.unsetf (ios :: floatfield);
    myfile5.precision(15);
    myfile5.open("shifted_to_zero.txt");

    for(int i=0;i<num;i++){
      myfile5 << shifted_to_zero[i][0] << " " << shifted_to_zero[i][1] << " " << shifted_to_zero[i][2] << " " << shifted_to_zero[i][3] << endl;
    }

    myfile5.close();

    double EVB[rows_common][cols+2];
    evb(shifted_to_autosub,EVB,num);

    myfile4.open("evb.txt"); 

    for(int i=0;i<num;i++){
      myfile4 << EVB[i][0] << " " << EVB[i][1] << " " << EVB[i][2] << " " << EVB[i][3] << endl;
    } 
    myfile4.close();
  }
  return 0;
}


void evb(double shifted_to_autosub[rows_common][cols+2],double EVB[rows_common][cols+2],int num){

  double delta_E[num];

//Tell Srini about this approximation!!!!!!!!!!! that some times delta_E_square is negative so you took modulus

  for(int i=0;i<num;i++){
    delta_E[i] = sqrt(pow((shifted_to_autosub[i][2] - shifted_to_autosub[i][0])*(shifted_to_autosub[i][2] - shifted_to_autosub[i][1]),2));

    //Ask Srini about using the unshifted values for individual topologies or the shifted ones
    //right now I am using the unshifted for the individual topologies
    if(shifted_to_autosub[i][3] != 0.0){
      EVB[i][0] = shifted_to_autosub[i][0];
      EVB[i][1] = shifted_to_autosub[i][1];
      EVB[i][2] = 0.5*(shifted_to_autosub[i][0] + shifted_to_autosub[i][1] - sqrt(pow((shifted_to_autosub[i][0] - shifted_to_autosub[i][1]),2) + 4*delta_E[i]));
      EVB[i][3] = shifted_to_autosub[i][3];
      
      cout << EVB[i][0] << " " << EVB[i][1] << " " << EVB[i][2] << " " <<  EVB[i][3] << endl;
    }
  }
}

void shifting_to_autosub(double shifted_to_zero[rows_common][cols+2],double shifted_to_autosub[rows_common][cols+2],double autosub_distance,int num){
  double shift_for_in = 0.0;
  double shift_for_trimer = 0.0;

  for(int i=0;i<num;i++){
    if(shifted_to_zero[i][3] == autosub_distance){
      shift_for_in = shifted_to_zero[i][0] - shifted_to_zero[i][1];
      shift_for_trimer = shifted_to_zero[i][2] - shifted_to_zero[i][1];
      cout << "Shift at autosub point for included topology is " << shift_for_in << endl;
      cout << "Shift at autosub point for trimer is " << shift_for_trimer << endl;
    }
  }

  for(int i=0;i<num;i++){
    shifted_to_autosub[i][0] = shifted_to_zero[i][0] - shift_for_in;
    shifted_to_autosub[i][1] = shifted_to_zero[i][1];
    shifted_to_autosub[i][2] = shifted_to_zero[i][2] - shift_for_trimer;
    shifted_to_autosub[i][3] = shifted_to_zero[i][3];
//    cout << shifted_to_autosub[i][0] << " " << shifted_to_autosub[i][1] << " " << shifted_to_autosub[i][2] << " " << shifted_to_autosub[i][3] << endl;
  }
}

void min(double common[rows_common][cols+2],double shifted_to_zero[rows_common][cols+2],int number_of_common_points){
  double min[3];
  min[0] = common[0][0];
  min[1] = common[0][1];
  min[2] = common[0][2];
  for(int i=0; i < number_of_common_points; i++){
    if(common[i][0] < min[0])
      min[0] = common[i][0];
    if(common[i][1] < min[1])
      min[1] = common[i][1];
    if(common[i][2] < min[2])
       min[2] = common[i][2];
  }
  cout.unsetf (ios :: floatfield);
  cout.precision(15);

  cout << min[0] << " " << min[1] << " " << min[2] << endl;

  for(int i=0; i < number_of_common_points; i++){
    if((common[i][0] != 0.0)&&(common[i][1] != 0.0)&&(common[i][2] != 0.0)){
      shifted_to_zero[i][0] = common[i][0] - min[0];
      shifted_to_zero[i][1] = common[i][1] - min[1];
      shifted_to_zero[i][2] = common[i][2] - min[2];
      shifted_to_zero[i][3] = common[i][3];
    }
    else
    cout << common[i][0] << " " << common[i][1] << " " << common[i][2] << " " << common[i][3] << i << endl;
  }
}

void common(double points1[rows1][cols],double points2[rows2][cols],double points3[rows3][cols],double common[rows_common][cols+2],int& number_of_common_points){
    double points[rows_common][cols];
    for(int i=0 ;i < rows1; i++){
        points[i][0] = points1[i][0];
        points[i][1] = points1[i][1];
    }

    for(int i=0; i< rows2; i++){
        points[i + rows1][0] = points2[i][0];
        points[i + rows1][1] = points2[i][1];
    }

    for(int i=0; i< rows3; i++){
        points[i + rows1 + rows2][0] = points3[i][0];
        points[i + rows1 + rows2][1] = points3[i][1];
   }

   ofstream myfile15;
   myfile15.unsetf (ios :: floatfield);
   myfile15.precision(15);
   myfile15.open("points.txt");

   for(int i=0; i< rows_common;i++){
     myfile15 << points[i][0] << " " << points[i][1] << " " << i << endl;
   }

   number_of_common_points = 0;
   ofstream myfile14;
   myfile14.unsetf (ios :: floatfield);
   myfile14.precision(15);
   myfile14.open("common.txt");

  for (int i = 0; i < rows_common; i++){
    for (int j = i + 1;j < rows_common; j++){
      for (int k = j + 1;k < rows_common; k++){
        if ((points[i][1] == points[j][1])&&(points[j][1] == points[k][1])&&(points[i][1] != 0)){
          common[number_of_common_points][0] = points[i][0];
          common[number_of_common_points][1] = points[j][0];
          common[number_of_common_points][2] = points[k][0];
          common[number_of_common_points][3] = points[i][1];
          myfile14 << common[number_of_common_points][0] << " " << common[number_of_common_points][1] << " " << common[number_of_common_points][2] << " " << common[number_of_common_points][3] << endl;
          number_of_common_points++;
        }
      }
    }
  }
}
