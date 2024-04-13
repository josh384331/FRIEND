#pragma once

#include <set>
#include "helper.hpp"

class Dataset
{
	
public:
    Dataset(){};
    ~Dataset(){};
    
    sarray mIndependentVarNames;
    sarray mDependentVarNames;

    void create_from_file(const string& filename);
    void create_from_data(dmatrix inputData, int numIndependentVars, vector<string> independentVarNames, vector<string> dependentVarNames);
    void print_data(int spacer);
    darray linear_interpolate(darray independentValues, int verbose);

private:
	//Private Variables
    string mFilename;
    int mNumIndependentVars;
    bool mHasSubset;

    vector<Dataset*> mSubset;  // Used if there is a subset of data
    darray mX; // Vector containing the independent variable values
    dmatrix mY; // vector containing the dependent variable values
    
};




