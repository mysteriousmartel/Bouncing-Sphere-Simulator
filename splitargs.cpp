// Copyright 2019 Jennifer Campbell janeeng@bu.edu
// Copyright 2019 David Henderson dth15@bu.edu

#include <fstream>
#include <iostream>

using namespace std;
// ofstream table;

int main(int argc, char** argv)
{
  // table.open("stdoutcpp.txt");

  for (int i = 1; i <= 4; ++i)
  {
    // table << *(argv + i) << "\n";
    cout << *(argv + i) << "\n";
  }

  // table.close();
  // table.open("stderrcpp.txt");

  for (int j = 5; j < argc; j++)
  {
    // table << *(argv + j) << "\n";
    cerr << *(argv + j) << "\n";
  }

  // table.close();

  return 0;

}
//i need to change the loop so that the file doesn't keep getting rewritten with each i