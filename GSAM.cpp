#include<vector>
#include<iostream>
#include<string>
#include<string.h>
#include <fstream>

using namespace std;
int GSAM[3000][3000];

int maximum(int a,int b,int c)
{
  int max = a;
  if(max<b) max=b;
  if(max<c) max = c;
  return max;
}
int calculatePointsDistance(string  string1, string string2 );

int main(int argc, char const *argv[]) {

  cout<<"hello"<<endl;

  ifstream file("./genomedata.csv");
  string line;
  vector<string> points;
  vector<string> sequences;
  string sequence="";

  while(getline(file,line))  //taking input
  {
      if(line[0] == '>')
      {
        string dataPoint ;
        for(int i=1;i<line.size();i++)
        {
          dataPoint.push_back(line[i]);
        }
        if(dataPoint.size()!=0)
        {
          points.push_back(dataPoint);
        }
        if(sequence.size()!=0)
        {
          sequences.push_back(sequence);
        }
        sequence = "";
      }
      else
      {
        for(int i=0;i<line.size();i++)
        {
          sequence.push_back(line[i]);
        }
      }

  }

  int pointsDistanceMatrix[points.size()][points.size()];
  for (size_t i = 0; i < points.size(); i++) { //initalize matrix to 0
    for (size_t j = 0; j < points.size(); j++) {
      pointsDistanceMatrix[i][j]=0;
    }
  }

  int counter = 0;
  for (size_t i = 0; i < points.size(); i++) { //initalize matrix to 0
    for (size_t j = i+1; j < points.size(); j++) {
      int dist = calculatePointsDistance(sequences[i],sequences[j]);
      //cout<<dist<<"sequence:" <<i<<" "<<j<<endl;
      pointsDistanceMatrix[i][j] = dist;
      pointsDistanceMatrix[j][i] = dist;
      counter+=1;

    }
    cout<<i<<endl;
  }

  ofstream opfile;
  opfile.open ("pointsDistanceMatrix.csv");
  for (size_t i = 0; i < points.size(); i++) {
    for (size_t j = 0; j < points.size(); j++) {
      opfile << pointsDistanceMatrix[i][j];
      if(j!=(points.size()-1))
      {
        opfile<<",";
      }
    }
    opfile<<"\n";
  }
  opfile.close();
  cout<<"end"<<endl;

  return 0;
}

int calculatePointsDistance(string  string1, string string2 )
{
    int gap= -2;
    int match = 1;
    int mismatch = -1;

    string sequence1;
    sequence1.push_back(' ');
    for (size_t i = 0; i < string1.size(); i++) {
      sequence1.push_back(string1[i]);
    }


    string sequence2;
    sequence2.push_back(' ');
    for (size_t i = 0; i < string2.size(); i++) {
      sequence2.push_back(string2[i]);
    }



    for (size_t i = 0; i < sequence1.size(); i++) {
      for (size_t j = 0; j < sequence2.size(); j++) {
        GSAM[i][j]=0;
      }
    }
    //cout << string1.size() << endl<< string2.size() << endl;

    int j,i;
    for (int i = 1; i < sequence1.size() ; i++) {
      GSAM[i][0] = GSAM[i-1][0] + gap;
    }
    for (int j = 1; j < sequence2.size(); j++) {
      GSAM[0][j] = GSAM[0][j-1] + gap;
    }

    for (int i = 1; i < sequence1.size(); i++) {
      for (int j = 1; j < sequence2.size(); j++) {
          int topChange = GSAM[i-1][j] + gap;
          int leftChange = GSAM[i][j -1] + gap;
          int diagonalChange;

          if( sequence1[i] == sequence2[j])
          {
            diagonalChange = GSAM[i-1][j-1] + match;
            //cout<<diagonalChange<<endl;
          }
          else
          {
            diagonalChange = GSAM[i-1][j-1] + mismatch;
          }

          GSAM[i][j] = maximum(topChange,leftChange,diagonalChange);


      }
    }

    return GSAM[sequence1.size()-1][sequence2.size()-1];

}
