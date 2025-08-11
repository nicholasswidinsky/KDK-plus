/**
 *
 * Code to sort coincidences in CAEN data.
 * P. Di Stefano, Queen's, 2023 09 25
 *
 * Paths are located in main() at end of file
 * Main function is sortData().  Read it to understand main features.
 *
 * To compile: 
 * g++ -o sortCAENdata sortCAENdata.cpp
 *
 * To execute:
 * ./sortCAENdata
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>	// For timing


using namespace std;


static double ps_ns = 0.001; // Conversion of ps to ns
static double nCh = 6; //16; // Number of channels
static int nDatCh = 3; // Number of output data items per channel (t, E, flag)



// Loads CSV data to a 2d array of strings, including headers.  Can't believe C++ doesn't have a library for this :(
vector<vector<string> > loadData(string fname)
{
    struct timeval begin, end;	// Timing
    gettimeofday(&begin, 0);	
    
	cout<<"Reading "<<fname<<"\n";
	
	vector<vector<string> > content;
	vector<string> row;
	string line, word;
 
	fstream file (fname, ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();
 
			stringstream str(line);
 
			while(getline(str, word, ';'))
				//cout<<word<<"\t";
				row.push_back(word);
			content.push_back(row);
		}
		//cout<<"\n";
	}
	else
		cout<<"Could not open the file\n";

	// Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);	//Timing
    long seconds = end.tv_sec - begin.tv_sec;
	cout << "\tReading time (s): " << seconds << endl;
	
 	
	return content;
}


/**
 * This function does the actual sorting into coins.
 *
 * First, creates the headers.  Each channel has 3 items (time, energy_AU, flag).
 * Then, creates a 2D array called "sortedData" with the same number of lines as the input data inData, and with all values set to the string 0.
 * Next, does the actual sorting.  For this, loops through inData line by line.  For each event, checks the following ones as long as they are in the coincident window.
 * Writes such events to the same line, in the appropriate channel column of sortedData.
 * Then, if clean is false, moves to the next event.  This means that coincident events will appear multiple times in sortedData.
 * If clean is true, then skips to the next uncoincident event.  Each event therefore only appears a single time, and  rows of sortedData that would be left at the default string 0 are removed.  However, the sum of all non-zero events in all the columns should be equal to the initial number of events.
 * 
 *
 * @return the sorted 2d array, with headers
 */
vector<vector<string> > sortData(vector<vector<string> > inData, double cwLo_ns, double cwHi_ns, bool clean, bool verbose)
{
    struct timeval begin, end;	// Timing
    gettimeofday(&begin, 0);	
    
	cout<<"Sorting"<<"\n";

	unsigned int n = inData.size();	// Includes header and data
	//unsigned int n = 10;	// Includes header and data
	
	vector<std::vector<string> > sortedData(n, std::vector<string>(nDatCh*nCh, "0"));	// Initialize sorted output to '0'; data will be E00_AU T00_ns Flag00 ... E15_AU T15_ns Flag15

	// Fill headers of 2d array sortedData
	for(unsigned int j = 0; j < nCh; j++)
	{
		string hdr = "T";
		hdr += to_string(j);
		hdr += "_ps";
		sortedData[0][0+j*nDatCh] = hdr;
		
		hdr = "E";
		hdr += to_string(j);
		hdr += "_AU";
		sortedData[0][1+j*nDatCh] = hdr;
		
		hdr = "Flag";
		hdr += to_string(j);
		sortedData[0][2+j*nDatCh] = hdr;
	}
    gettimeofday(&end, 0);	//Timing
    long seconds = end.tv_sec - begin.tv_sec;
	cout << "\tInitialization time (s): " << seconds << endl;
	

	
	// The main part
	// Loop through rows of data and find coins
	unsigned int nCoin = 0;	// Number of coincidences that have been found (of any, including single, multiplicity).  This is used to fill the output file.
	
	for (unsigned int i = 1; i < n; i++)	// Start at 1 to skip header
	{	
		nCoin ++;
		
		if (((i<=30) & ((i%10)==0)) || ((i<=300) & ((i%100)==0)) || ((i<=3000) & ((i%1000)==0)) || ((i<=30000) & ((i%10000)==0)) || ((i<=300000) & ((i%100000)==0)) || ((i%300000)==0)) {
			gettimeofday(&end, 0);	//Timing
			long seconds = end.tv_sec - begin.tv_sec;
			cout << "\t\tEvent " << i << ", sorting time (s): " << seconds << endl;
			
		}
		
		// Search forward
		unsigned int j = i + 1;
		double dt_ns = 0;
		
		int ch1 = stoi(inData[i][1]);
		long long int t1_ps = stoll(inData[i][2]);
		string E1_AU = inData[i][3];
		string flag1 = inData[i][5];
		
		sortedData[nCoin][3*ch1+0] = to_string(t1_ps);	// Fill time
		sortedData[nCoin][3*ch1+1] = E1_AU;	// Fill energy
		sortedData[nCoin][3*ch1+2] = flag1;	// Fill string
		
		if (verbose) {
			cout<<"i="<<i<<"\tCh="<<ch1<<"\tt (ps)="<<t1_ps<<endl;
		}
		
		while ((j < n) && (dt_ns >= cwLo_ns) && (dt_ns <= cwHi_ns)) {
			int chN = stoi(inData[j][1]);
			long long int tN_ps = stoll(inData[j][2]);
			string EN_AU = inData[j][3];
			string flagN = inData[j][5];
			
			dt_ns = (tN_ps - t1_ps) * ps_ns;
			
			if ((dt_ns >= cwLo_ns) && (dt_ns <= cwHi_ns)) {
				sortedData[nCoin][3*chN+0] = to_string(tN_ps);	// Fill time
				sortedData[nCoin][3*chN+1] = EN_AU;	// Fill energy
				sortedData[nCoin][3*chN+2] = flagN;	// Fill string
				
				if (clean) {
					i = j; // Increment i to avoid events being counted multiple times
				}
			}
			
			if (verbose) {
				cout<<"ooooooooooooooooooooo "<<"j="<<j<<"\tCh="<<chN<<"\tt (ps)="<<tN_ps<<"\tdt (ns)="<<dt_ns<<endl;
			}
			
			j++;
			/* 
			I am putting all the variables that I am deleting under here. If things catastrophically break I know where they all are. 
			*/




		}
	}
    gettimeofday(&end, 0);	//Timing
    seconds = end.tv_sec - begin.tv_sec;
	cout << "\tSorting time (s): " << seconds << endl;
	


	// Remove empty lines.
	if (clean) {
		sortedData.erase(sortedData.begin()+nCoin+1, sortedData.end());		
	}

    gettimeofday(&end, 0);	//Timing
    seconds = end.tv_sec - begin.tv_sec;
	cout << "\tCleaning time (s): " << seconds << endl;

	
	return sortedData;
}


// Write data to tab-delimited csv file.  Can't believe C++ doesn't have a library for this :(
void writeData(vector<vector<string> > inData, string fname)
{
    struct timeval begin, end;	// Timing
    gettimeofday(&begin, 0);	
    
	cout<<"Writing to "<<fname<<endl;


	std::ofstream myfile;
	myfile.open(fname);

	unsigned int n = inData.size();
	for(unsigned int i = 0; i < n; i++)
	{
		for(unsigned int j = 0; j < nDatCh*nCh; j++)
		{
			myfile << inData[i][j];
			if (j < nDatCh*nCh-1) {
				myfile << "\t";
			}
			else if (i < n-1) {
				myfile << "\n";
			}
		}
	}
	
	myfile.close();

    gettimeofday(&end, 0);	//Timing
    long seconds = end.tv_sec - begin.tv_sec;
	cout << "\tWriting time (s): " << seconds << endl;
	
	
}

/**
 *
 *
 * Optional arguments:
 *
 * inName
 * outName
 *
 * cwLo_ns
 * cwHi_ns
 *
 *
 * clean
 * verbose
 */
int main(int argc, char *argv[])
{
	bool clean = true;

	bool verbose = false;
	
	double cwLo_ns = -300;
	double cwHi_ns = 300;

	string inName = "/home/nick/KDK+/LLT_3ml_1M/RAW/SDataR_LLT_3ml_1M.CSV";
	string outName = "/home/nick/KDK+/LLT_3ml_1M/RAW/new_SDataR_LLT_3ml_1M_coinSorted.csv";
	//string inName = "/Users/distefano/Documents/Data/TestData/CAEN_Test/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1200_200lsb_15h_raw/RAW/SDataR_2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1200_200lsb_15h_raw.CSV";
	//string outName = "/Users/distefano/Documents/Data/TestData/CAEN_Test/2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1200_200lsb_15h_raw/RAW/SDataR_2SiPM_coincidence_80d_50d_NoCoincBoard_240nsCoMPASS_CoincidenceWindow_1200_200lsb_15h_raw_COINSORTED.CSV";

	if (argc>=3) {
		inName = argv[1];
		outName = argv[2];
	}
	if (argc>=5) {
		cwLo_ns = stod(argv[3]);
		cwHi_ns = stod(argv[4]);
	}
	if (argc>=7) {
		clean = stoi(argv[5]);
		verbose = stoi(argv[6]);
	}
	
	
	vector<vector<string> > input = loadData(inName);
	
	vector<vector<string> > sorted = sortData(input, cwLo_ns, cwHi_ns, clean, verbose);
	input.clear();
	writeData(sorted, outName);
	
	cout<<"That's all folks!"<<endl;
	return 0;
}
 