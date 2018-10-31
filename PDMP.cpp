// A stochastic and dynamical view of pluripotency in mouse embryonic stem cells
// Published in PLOS Computational Biology

// Authors: 
// Yen Ting Lin (yentingl@lanl.gov, the one who coded this up and should be responsible...)
// T-6 and T-CNLS, Los Alamos National Laboratory, Los Alamos, NM 87544, USA
// Peter G. Hufton (co-developer)
// School of Physics and Astronomy, The University of Manchester, M13 9PL, UK
// Esther J. Lee, Department of Bioengineering, Rice University, Houston, TX 77005, USA
// Davit A. Potoyan, Department of Chemistry, Iowa State University, Ames, IA 50011, USA

// Instruction, compile the c++ source
// > g++ -o a.out PDMP.cpp
// Then, the executable takes 6 arguments: (LIF, CH, PD) before induction, (LIF, CH, PD) after induction
// For example, if initial LIF=0, CH=0, PD=1, and after induction LIF=1, CH=1, PD=1, execute
// > ./a.out 0 0 1 1 1 1
// The output file will be the marginal probability distributions with format specified below (**)

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <iomanip>

using namespace std;

struct parameter
{
	int numberGenes;	
	// total number of the genes, call it N
	
	int *geneType;			
	// type of the gene. 0=only regulated by activators, 1=only regulated by repressors, and 2=both
	
	double *basalTranscriptionRate;			
	// Basal rate of type-i gene is basalRate[i].
	
	double *gamma;		
	// The degradation rate of the ith type of TF is gamma[i].	
	
	double **diffTranscriptionRate;		
	// The production rate of the TF, NxNx2 matrix. When ith type of the TF binds to the jth gene, the transcription rate of the gene is increased by beta[i][j].
	
	double **koff; 		
	// The dissociation rate of bound (type-i) TF to (type-j) gene is equal to koff[i][j].
	
	double **kon; 		
	// The (per-capita) binding rate of free (type-i) TF to (type-j) gene is equal to kon[i][j].
	
		
	int **Epos;
	// The positively regulating edge indices
	int dimEposx, dimEposy;
	// size of Epos
	
	int **Eneg;
	// The negative regulating edge indices
	int dimEnegx, dimEnegy;
	// size of Epos
	
	int **Upos;
	// The upstream positive regulatory index
	int dimUposx, dimUposy;
	// size of Upos;
	
	int **Uneg;
	// The upstream negative regulatory index
	int dimUnegx, dimUnegy;
	// size of Uneg;
	
	int **Dpos;
	// The downstream positive regulatory index
	int dimDposx, dimDposy;
	// size of Dpos;
	
	int **Dneg;
	// The downstream negative regulatory index
	int dimDnegx, dimDnegy;
	// size of Dneg;
	
	
	int *nDpos;
	// the number of downstream, positively regulating gene
	int *nDneg;
	// the number of downstream, negatively regulating gene
	int *nUpos;
	// the number of upstream, positively regulating gene
	int *nUneg;
	// the number of upstream, negatively regulating gene
	
	int Nmax;
	// maximum binding sites
	
	double LIF, CH, PD, uniformTranscriptionRate,uniformGamma,uniformKoff,uniformKon,uniformLeak,uniformBasal;
	
};

struct state
{
	
	double currentTime;
	// the current time
	
	double nextTime;
	// time at which next switching event occurs
	
	int nextIndex;
	// index of the next switching event
	
	double *nextSwitchingTime;
	// time at which next switching event of gene i occurs
	
	int *nextSwitchingIndex;
	// index of the next switching event of gene i (labelled by the edge map)
	
	double *TFDensity;
	// the density of the TF
	
	int *geneticState;
	// the state of the gene. 0: unbound, -1: bound to a repressor, and 1: bound to a activator
	
	int *currentlyBoundTFIndex;
	// the type of the TF bound to the gene; when it is not bound to a TF, the value is -1
	
	double *currentTranscriptionRate;
	// the transcription rate of the gene at the moment
	
	int **currentlyBoundTFArray;
	// used for multiple binding sites.  (i=1:15) and (j=1:2)
	// currentlyBoundTFArray[i][0]: the number of bound activators for gene i
	// currentlyBoundTFArray[i][1]: the number of bound repressors for gene i
	
	
	int **nextReactionDir;
	// used for multiple binding sites.  (i=1:15)
	// it deposits the direction of next move, format similar to currentlyBoundTFArray
	
};

void initiate_parameters(parameter * par);
void initiate_state(state * var, parameter par);
void recycle_state(state * var, parameter par);
double evolve_until_T(state *sta, parameter par, double tend);
double rndexclusive();
double exponential(double rate);
double waitingTimeSwitchON(double iniRate, double saturatingRate, double expRate);
double waitingTimeSwitchOFF(double rate);
double evolveTFDynamics(double x0, double saturation, double expRate, double dt);
void regenerateWaitingTimes(state *var, parameter par, int targetIndex);
void testVisualizingParameterSet(parameter par);
void testVisualizingVariable(state variable, parameter par);
void testWaitingTimeDistributions();

int main(int argc, char *argv[])
{
    srand (time(NULL));
	
	/**************************/
	/*initiation of parameters*/
	/**************************/
	int ensN = 1E5; // ensemble number;
	int Ngrid = 100; // # of histogram bins
	
	parameter par;
	
	
	// initial test: quenched LIF, CH, PD, and uniform rates
	par.Nmax = 2;
	par.LIF = 1;
	par.CH = 1;
	par.PD = 1;
	par.uniformGamma = 1.0;         // degradation rate gamma in paper
	par.uniformKon = 16.0;          // k_off in paper
	par.uniformKoff = 1.5;          // k_on in paper
	par.uniformLeak = 0.0;          // just let it be 0; in paper we did not present this
	par.uniformBasal = 0.01;        // alpha_m in paper
    
	par.uniformTranscriptionRate = 1;
	
	double LIF = atof(argv[1]);
	double CH = atof(argv[2]);
	double PD = atof(argv[3]);
	double LIF2 = atof(argv[4]);
	double CH2 = atof(argv[5]);
	double PD2 = atof(argv[6]);
	
	par.LIF = LIF;
	par.CH = CH;
	par.PD = PD;
	
	initiate_parameters(&par);
	
	state * var = new state [ensN];

	for (int i=0;i<ensN;i++)
		initiate_state(var+i, par);

	int G;
	double * histogram = new double [Ngrid];
	
	double dd = 1.0/Ngrid; // grid size, for computing proper probability measure
	
	ofstream output;
	stringstream filename;
	filename << "PDMPDyanmics" << LIF << CH << PD << "to" << LIF2 << CH2 << PD2 << ".txt";
	output.open(filename.str().c_str());
		
	for (double t= 0.0;t<20; t+=0.05)
	{
		
        // thermalization 
        
		for (int index=0;index <ensN;index++)
		{
			var[index].TFDensity[12]=LIF;
			var[index].TFDensity[13]=CH;
			var[index].TFDensity[14]=PD;
		}
		
		cout << "Evolving t=" << t << ", Tend=" << 60 << endl;
		
		for (int index = 0; index < ensN; index++)
		{
			evolve_until_T(var+index, par, t);
		}
		
		// here we output the config
		// each row is a "snapshot" for a time point
		// First column: time
		// 2 to end: record the distributions

		if (t>19)
		{
            
            // We only record this from t=19 to t=20 of the uninduced condition
            // t<19 discarded -- thermalization
            
			output << var[0].currentTime << "\t";
			
			for (int index2 = 0; index2 < 12; index2++)
			{
								
				for (int iii = 0;iii<Ngrid;iii++)
					histogram[iii] = 0.0;
					
					
				for (int index = 0; index < ensN; index++)
				{
					histogram[(int)floor(var[index].TFDensity[index2] / dd )]++;
				}	
				
				for (int kkk = 0;kkk<Ngrid;kkk++)
				{
					output << histogram[kkk]/ensN/dd << "\t";
				}
				
			}
			
			output << endl;
		}
		
	}

    // Induction , change external environment to new LIF, CH, PD
    
    par.LIF = LIF2;
    par.CH = CH2;
    par.PD = PD2;

    initiate_parameters(&par);
    
    for (int index=0;index <ensN;index++)
    {
        recycle_state(var+index, par);
    }
        
	for (double t= 20.0;t<60; t+=0.05)
	{
		
		cout << "Evolving t=" << t << ", Tend=" << 60 << endl;
		
		for (int index = 0; index < ensN; index++)
		{
			evolve_until_T(var+index, par, t);
		}
		
		// (**) here we output the config
		// each row is a "snapshot" for a time point
		// First column: time
		// 2 to end: the marginal distributions of each TF from density from 0 to 1, binned into Ngrid parts
        // Format is [P_{TF1 density in first bin}, P_{TF_1 density in second bin}...P_{TF_1 density in last bin}, P_{TF2 density in first bin}, P_{TF_2 density in second bin}...P_{TF_2 density in last bin},..., P_{TF_12 density in first bin}...P_{TF_12 density in last bin}]

        output << var[0].currentTime << "\t";
        
        for (int index2 = 0; index2 < 12; index2++)
        {
                            
            for (int iii = 0;iii<Ngrid;iii++)
                histogram[iii] = 0.0;
                
                
            for (int index = 0; index < ensN; index++)
            {
                histogram[(int)floor(var[index].TFDensity[index2] / dd )]++;
            }	
            
            for (int kkk = 0;kkk<Ngrid;kkk++)
            {
                output << histogram[kkk]/ensN/dd << "\t";
            }
            
        }
        
        output << endl;

    }
			
			
			
	return 0;

}

double evolveTFDynamics(double x0, double saturation, double expRate, double dt)
{
	double xf = saturation - (saturation - x0) * exp(-expRate*dt);
	return xf;
}

void testWaitingTimeDistributions()
{
	
	int N = 100;
	int ensN = 1E7;
	
	double tmax = 10.0;
	
	double ri = 0.0;
	double rs = 1.0;
	double gamma = 1.0;
	
	double * tgrid = new double [N];
	double * survival = new double [N];
	
	for (int i=0;i<N;i++)
	{
		tgrid[i] = (double)i*tmax/N;
		survival[i] = 0;
	}
	
	double dt ;
	
	int j = 0;
	
	for (int i=0;i<ensN;i++)
	{
		dt = waitingTimeSwitchON(ri, rs, gamma);
		//cout << dt << endl;
		
		j = 0;
		
		bool flag = true;
		
		while (flag)
		{
			
			survival[j] ++;
			j ++;
			flag = (tgrid[j]<dt)&&(j<N);
	
		}
	}
	
	ofstream output;
	output.open("TestWaitingTimes.txt");
	
	for (int i=0;i<N;i++)
	{
		output << tgrid[i] << "\t" << survival[i] / ensN << endl;
	}
	
	
	output.close();
}

void initiate_state(state * var, parameter par)
{
	
	(*var).currentTime = 0;
	(*var).nextTime = 0;
	(*var).nextIndex = -1;
	(*var).nextSwitchingTime = new double [par.numberGenes];
	(*var).nextSwitchingIndex = new int [par.numberGenes];
	(*var).TFDensity = new double [par.numberGenes];
	(*var).geneticState = new int [par.numberGenes];
	(*var).currentlyBoundTFIndex = new int [par.numberGenes];
	(*var).currentTranscriptionRate = new double [par.numberGenes];
	
	for (int i=0;i<par.numberGenes;i++)
	{
		(*var).nextSwitchingTime[i] = 1E20;
		(*var).nextSwitchingIndex[i] = -1;
		(*var).TFDensity[i] = 0.5;
		(*var).geneticState[i] = 0;
		(*var).currentlyBoundTFIndex[i] = -1;
		(*var).currentTranscriptionRate[i] = par.basalTranscriptionRate[i];
	}
	
	// Here we assign the constant LIF, CH, PD densities
	(*var).TFDensity[12] = par.LIF;		// LIF
	(*var).TFDensity[13] = par.CH;		// CH
	(*var).TFDensity[14] = par.PD;		// PD
	
	// now propose the next switching for each gene
	
	(*var).currentlyBoundTFArray = new int * [15];
	(*var).nextReactionDir = new int * [15];
	
	for (int i=0;i<15;i++)
	{
		(*var).currentlyBoundTFArray[i] = new int [2];
		(*var).nextReactionDir[i] = new int [2];
		
		(*var).currentlyBoundTFArray[i][0] = 0;
		(*var).currentlyBoundTFArray[i][1] = 0;	
		
		(*var).nextReactionDir[i][0] = 0;
		(*var).nextReactionDir[i][1] = 0;	
	}	
	
	

	for (int i=0;i<par.numberGenes;i++)
		regenerateWaitingTimes(var, par, i);
			
}

void recycle_state(state * var, parameter par)
{
	
	//(*var).currentTime = 0;
	//(*var).nextTime = 0;
	(*var).nextIndex = -1;
	
	for (int i=0;i<par.numberGenes;i++)
	{
		(*var).nextSwitchingTime[i] = 1E20;
		(*var).nextSwitchingIndex[i] = -1;
		(*var).currentTranscriptionRate[i] = par.basalTranscriptionRate[i];
	}
	
	// Here we assign the constant LIF, CH, PD densities
	(*var).TFDensity[12] = par.LIF;		// LIF
	(*var).TFDensity[13] = par.CH;		// CH
	(*var).TFDensity[14] = par.PD;		// PD
	
	// now propose the next switching for each gene
	
	for (int i=0;i<par.numberGenes;i++)
		regenerateWaitingTimes(var, par, i);
			
}

double waitingTimeSwitchON(double iniRate, double saturatingRate, double expRate)
{
	// this function generates the waiting times to switch on the gene
	// assuming a time-dependent rate having the form
	// r(t) = saturatingRate - (saturatingRate - iniRate) exp( - expRate * t).
	
	double dt = 0;
	
	if ( iniRate > saturatingRate)
	{
		// in this case, the rate decays and converges to saturatingRate
		// We generate two random waiting times: one for constant saturatingRate/
		// And another for exponentially decaying rate (in the note: Eq. 8)
		
		double dt1, dt2;
		if (saturatingRate>0)
		{
			dt1 = exponential(saturatingRate);
			
		}else{
			
			dt1 = 1E10;  // never gonna happen
			
		}
		
		double u = rndexclusive();
		
		if (u > exp(- (iniRate-saturatingRate)/expRate  ))
		{
			
			dt2 = -log(expRate/(iniRate-saturatingRate)*log(u) +1) / expRate;
			
		}else{
			
			dt2 = 1E10;
			
		}
		
		dt = min(dt1, dt2);
		
		
	}
	else if (iniRate == saturatingRate)
	{
		
		dt = exponential(iniRate);
		
	}
	else
	{
		// in this case, the rate increases and converges to saturatingRate 
		// We have to solve the inversion sampling problem numerically (i.e., Eq. 10)
	
        // change to your favorite random number generator!!
    	double u = ((double) rand() / (RAND_MAX));
		
		double LB = -log(u)/saturatingRate;
		double UB = LB + (saturatingRate-iniRate) / saturatingRate / expRate;
		double C = UB;
	
		double initialGuess = 1E10;
		double secondGuess = LB/2+UB/2;
		
		while (abs(initialGuess-secondGuess) > 1E-8) // tolerance 1E-6
		{
			initialGuess = secondGuess;
			// Newton's method to solve t + exp(-expRate*t) - C = 0;
			secondGuess = initialGuess -  (initialGuess + (saturatingRate - iniRate)/expRate/saturatingRate*exp(-expRate*initialGuess) - C) / (1 - (saturatingRate-iniRate)/saturatingRate*(-expRate*initialGuess));
			
		}

		dt = initialGuess / 2 + secondGuess / 2;
		

		if ((dt > UB)||(dt < LB))
			cout << "warning: the waiting time is not correct!!";
	}
	
	return dt;
	
	
}

double evolve_until_T(state *sta, parameter par, double tend)
{
	
	int geneIndex, geneIndexRoot, TFIndex;
	double currentRate, saturatingRate, expRate;
	bool switchingToON;
	
	int candidate = -1;
	
	//cout << "Current Time=" << (*sta).currentTime << endl;
	
	while ((*sta).currentTime < tend)
	{
		double proposedTime=1E10;
		candidate = -1;
		
		for (geneIndex = 0;geneIndex <par.numberGenes;geneIndex ++)
		{
			if ((*sta).nextSwitchingTime[geneIndex] < proposedTime)
			{
				proposedTime = (*sta).nextSwitchingTime[geneIndex];
				candidate = geneIndex;
			}
		}
		
		
		if (proposedTime < tend)
		{
			
			// Evolve the TFDensity

			
			for (int i=0;i<par.numberGenes;i++)
			{
								
				(*sta).TFDensity[i] = evolveTFDynamics((*sta).TFDensity[i],
													   (*sta).currentTranscriptionRate[i] / par.gamma[i],
													   par.gamma[i],									   
													   proposedTime - (*sta).currentTime);
				
			}		

			(*sta).currentTime = proposedTime;
			geneIndex = candidate;
			// update bound array
			
			(*sta).currentlyBoundTFArray[geneIndex][0] += (*sta).nextReactionDir[geneIndex][0];
			(*sta).currentlyBoundTFArray[geneIndex][1] += (*sta).nextReactionDir[geneIndex][1];
					
			// update transcription rate
			if (par.geneType[geneIndex]==0)
			{
				if ((*sta).currentlyBoundTFArray[geneIndex][0]==par.Nmax)
				{
					(*sta).currentTranscriptionRate[geneIndex] = par.uniformTranscriptionRate;
				}else{
					(*sta).currentTranscriptionRate[geneIndex] = par.uniformLeak;
				}
			}else if (par.geneType[geneIndex]==1){
				if ((*sta).currentlyBoundTFArray[geneIndex][1]==par.Nmax)
				{
					(*sta).currentTranscriptionRate[geneIndex] = par.uniformLeak;
				}else{
					(*sta).currentTranscriptionRate[geneIndex] = par.uniformTranscriptionRate;
				}
			}else{
				if ((*sta).currentlyBoundTFArray[geneIndex][0]==par.Nmax)
				{
					(*sta).currentTranscriptionRate[geneIndex] = par.uniformTranscriptionRate;
				}else if ((*sta).currentlyBoundTFArray[geneIndex][1]==par.Nmax){
					(*sta).currentTranscriptionRate[geneIndex] = par.uniformLeak;
				}else{
					(*sta).currentTranscriptionRate[geneIndex] = par.uniformBasal;
				}
			}
		
			
			// update the random waiting times;
						
			geneIndexRoot = geneIndex;
						
			regenerateWaitingTimes(sta, par, geneIndexRoot);
			
		
			// As the gene "geneIndex" changes, the downstream waiting times must be re-drawn
			for (int kk = 0;kk < par.nDpos[geneIndexRoot]; kk++)
				regenerateWaitingTimes(sta, par, par.Epos[1][par.Dpos[kk][geneIndexRoot]]);
			
			for (int kk = 0;kk < par.nDneg[geneIndexRoot]; kk++)
				regenerateWaitingTimes(sta, par, par.Eneg[1][par.Dneg[kk][geneIndexRoot]]);
			
			
		}else{
			
			
			// No switching event occur before tend
			
			(*sta).nextTime = tend;

			// Evolve the TFDensity

			for (int i=0;i<par.numberGenes;i++)
			{
				
				(*sta).TFDensity[i] = evolveTFDynamics((*sta).TFDensity[i],
													   (*sta).currentTranscriptionRate[i] / par.gamma[i],
													   par.gamma[i],									   
													   (*sta).nextTime - (*sta).currentTime);
				
			}		
			
			(*sta).currentTime = (*sta).nextTime;			
				
		}
		
		
	}
	
}

void regenerateWaitingTimes(state *sta, parameter par, int geneIndex)
{
			
			
	// this function regenerates all the waiting times for possible binding reactions on geneIndex
	int TFIndex;
	double currentRate, saturatingRate, dt;
	(*sta).nextSwitchingTime[geneIndex] = 1E20;
	int totalBoundTF =(*sta).currentlyBoundTFArray[geneIndex][0]+(*sta).currentlyBoundTFArray[geneIndex][1];
	
	
	if (totalBoundTF==0)
	{
		
		
		// All sites empty, can only bind to a transcription factor
		
		for (int ll=0;ll<par.nUpos[geneIndex]; ll++)
		{
	 

			TFIndex = par.Epos[0][par.Upos[ll][geneIndex]];
			
			currentRate = par.Nmax * par.uniformKon * (*sta).TFDensity[TFIndex];
			saturatingRate = par.Nmax * par.uniformKon * (*sta).currentTranscriptionRate[TFIndex] / par.gamma[TFIndex];

			dt = waitingTimeSwitchON(currentRate, saturatingRate, par.gamma[TFIndex]);


			if ((*sta).currentTime + dt < (*sta).nextSwitchingTime[geneIndex])
			{

				(*sta).nextSwitchingTime[geneIndex] = (*sta).currentTime + dt;
				(*sta).nextReactionDir[geneIndex][0] = +1;
				(*sta).nextReactionDir[geneIndex][1] = 0;

			}		
			
		}		

		
		for (int ll=0; ll < par.nUneg[geneIndex]; ll++)
		{
			
			TFIndex = par.Eneg[0][par.Uneg[ll][geneIndex]];
			
			currentRate = par.Nmax * par.uniformKon * (*sta).TFDensity[TFIndex];
			saturatingRate = par.Nmax * par.uniformKon * (*sta).currentTranscriptionRate[TFIndex] / par.gamma[TFIndex];

			dt = waitingTimeSwitchON(currentRate, saturatingRate, par.gamma[TFIndex]);
			


			if ((*sta).currentTime + dt < (*sta).nextSwitchingTime[geneIndex])
			{

				(*sta).nextSwitchingTime[geneIndex] = (*sta).currentTime + dt;
				(*sta).nextReactionDir[geneIndex][0] = 0;
				(*sta).nextReactionDir[geneIndex][1] = +1;

			}		
			
		}					
		
		
	
	}else if (totalBoundTF == par.Nmax){ 
		
		
		// All sites occupied, can only unbind
		
		(*sta).nextSwitchingTime[geneIndex] = (*sta).currentTime + waitingTimeSwitchOFF(par.Nmax*par.uniformKoff);
		
        // change to your favorite random number generator!!
    	double u = ((double) rand() / (RAND_MAX))*par.Nmax;
		
		
		if (u<(*sta).currentlyBoundTFArray[geneIndex][0])
		{
			(*sta).nextReactionDir[geneIndex][0] = -1;
			(*sta).nextReactionDir[geneIndex][1] = 0;
		}else{
			(*sta).nextReactionDir[geneIndex][0] = 0;
			(*sta).nextReactionDir[geneIndex][1] = -1;			
		}
		
		
	}else{
		// in between, can do both
 		
		for (int ll=0;ll<par.nUpos[geneIndex]; ll++)
		{
	 

			TFIndex = par.Epos[0][par.Upos[ll][geneIndex]];
			
			currentRate = (par.Nmax-totalBoundTF) * par.uniformKon * (*sta).TFDensity[TFIndex];
			saturatingRate = (par.Nmax-totalBoundTF) * par.uniformKon * (*sta).currentTranscriptionRate[TFIndex] / par.gamma[TFIndex];

			dt = waitingTimeSwitchON(currentRate, saturatingRate, par.gamma[TFIndex]);


			if ((*sta).currentTime + dt < (*sta).nextSwitchingTime[geneIndex])
			{

				(*sta).nextSwitchingTime[geneIndex] = (*sta).currentTime + dt;
				(*sta).nextReactionDir[geneIndex][0] = +1;
				(*sta).nextReactionDir[geneIndex][1] = 0;

			}		
			
		}		

		
		for (int ll=0; ll < par.nUneg[geneIndex]; ll++)
		{
			
			TFIndex = par.Eneg[0][par.Uneg[ll][geneIndex]];
			
			currentRate = (par.Nmax-totalBoundTF) * par.uniformKon * (*sta).TFDensity[TFIndex];
			saturatingRate = (par.Nmax-totalBoundTF) * par.uniformKon * (*sta).currentTranscriptionRate[TFIndex] / par.gamma[TFIndex];

			dt = waitingTimeSwitchON(currentRate, saturatingRate, par.gamma[TFIndex]);
			


			if ((*sta).currentTime + dt < (*sta).nextSwitchingTime[geneIndex])
			{

				(*sta).nextSwitchingTime[geneIndex] = (*sta).currentTime + dt;
				(*sta).nextReactionDir[geneIndex][0] = 0;
				(*sta).nextReactionDir[geneIndex][1] = +1;

			}		
			
		}		

		dt = (*sta).currentTime + waitingTimeSwitchOFF(totalBoundTF*par.uniformKoff);
		
		if (dt < (*sta).nextSwitchingTime[geneIndex] )
		{
	
            // change to your favorite random number generator!!
            double u = ((double) rand() / (RAND_MAX))*totalBoundTF;
			
			if (u<(*sta).currentlyBoundTFArray[geneIndex][0])
			{
				(*sta).nextReactionDir[geneIndex][0] = -1;
				(*sta).nextReactionDir[geneIndex][1] = 0;
			}else{
				(*sta).nextReactionDir[geneIndex][0] = 0;
				(*sta).nextReactionDir[geneIndex][1] = -1;			
			}
			
		}		
		

		
	}
	
	
}

double waitingTimeSwitchOFF(double rate)
{
	// this function generate waiting times for dissociation events
	if (rate==0)
	{
		return 1E20;
	}else{	
		double dt = exponential(rate);
		return dt;
	}
}

double exponential(double rate)
{
	double out = -1.0/(rate) * log(rndexclusive());
	return out;
}

double rndexclusive()
{
    
    // change to your favorite random number generator!!
    
	double i = 0;
    
    while ((i==0)||(i==1))
	{
		i = ((double) rand() / (RAND_MAX));
	}
    
    return i;
}  

void initiate_parameters(parameter * par)
{
	
	// degradation rate is set to be 1.0 (suitable timescale and population renormalised by some large population scale)
	
	// (*par).LIF = 1.0;
	// (*par).CH = 1.0;
	// (*par).PD = 1.0;
	
	(*par).numberGenes = 15;
	
	(*par).geneType = new int [(*par).numberGenes];
		
	(*par).basalTranscriptionRate = new double [(*par).numberGenes];
	(*par).gamma    			  = new double [(*par).numberGenes];
	
	
	(*par).diffTranscriptionRate = new double * [(*par).numberGenes];
	(*par).koff 				 = new double * [(*par).numberGenes];
	(*par).kon  				 = new double * [(*par).numberGenes];
		
	
	for (int i=0;i<(*par).numberGenes;i++)
	{
		(*par).diffTranscriptionRate[i] = new double [(*par).numberGenes];
		(*par).koff[i] 				    = new double [(*par).numberGenes];
		(*par).kon[i] 				    = new double [(*par).numberGenes];
		
		for (int j=0;j<(*par).numberGenes;j++)
		{
			(*par).diffTranscriptionRate[i][j]  = 0.0;
			(*par).koff[i][j] 				    = 0.0;
			(*par).kon[i][j] 				    = 0.0;
		}
		
	}
	
		
	// Setting up the positive edge index
	(*par).Epos = new int * [2];
	(*par).Epos[0] = new int [18];
	(*par).Epos[1] = new int [18];
	(*par).dimEposx = 2;
	(*par).dimEposy = 18;
	
	(*par).Epos[0][0] = 0;		
	(*par).Epos[0][1] = 0;
	(*par).Epos[0][2] = 0;
	(*par).Epos[0][3] = 3;
	(*par).Epos[0][4] = 4;
	(*par).Epos[0][5] = 4;
	(*par).Epos[0][6] = 5;
	(*par).Epos[0][7] = 6;
	(*par).Epos[0][8] = 6;
	(*par).Epos[0][9] = 7;
	(*par).Epos[0][10] = 8;
	(*par).Epos[0][11] = 8;	
	(*par).Epos[0][12] = 9;	
	(*par).Epos[0][13] = 9;	
	(*par).Epos[0][14] = 10;	
	(*par).Epos[0][15] = 11;	
	(*par).Epos[0][16] = 11;	
	(*par).Epos[0][17] = 12;	

	(*par).Epos[1][0] = 3;		
	(*par).Epos[1][1] = 4;
	(*par).Epos[1][2] = 6;
	(*par).Epos[1][3] = 6;
	(*par).Epos[1][4] = 5;
	(*par).Epos[1][5] = 9;
	(*par).Epos[1][6] = 4;
	(*par).Epos[1][7] = 4;
	(*par).Epos[1][8] = 11;
	(*par).Epos[1][9] = 8;
	(*par).Epos[1][10] = 5;
	(*par).Epos[1][11] = 10;	
	(*par).Epos[1][12] = 10;	
	(*par).Epos[1][13] = 11;	
	(*par).Epos[1][14] = 7;	
	(*par).Epos[1][15] = 7;	
	(*par).Epos[1][16] = 8;	
	(*par).Epos[1][17] = 0;		
	
	
	
	// Setting up the eegative dge index
	(*par).Eneg = new int * [2];
	(*par).Eneg[0] = new int [8];
	(*par).Eneg[1] = new int [8];
	(*par).dimEnegx = 2;
	(*par).dimEnegy = 8;
	
	(*par).Eneg[0][0] = 1;		
	(*par).Eneg[0][1] = 1;
	(*par).Eneg[0][2] = 2;
	(*par).Eneg[0][3] = 2;
	(*par).Eneg[0][4] = 5;
	(*par).Eneg[0][5] = 7;
	(*par).Eneg[0][6] = 13;
	(*par).Eneg[0][7] = 14;

	(*par).Eneg[1][0] = 4;		
	(*par).Eneg[1][1] = 5;
	(*par).Eneg[1][2] = 1;
	(*par).Eneg[1][3] = 8;
	(*par).Eneg[1][4] = 7;
	(*par).Eneg[1][5] = 4;
	(*par).Eneg[1][6] = 1;
	(*par).Eneg[1][7] = 2;
	
	
	// Setting up the Dpos index
	(*par).Dpos = new int * [3];
	(*par).Dpos[0] = new int [15];
	(*par).Dpos[1] = new int [15];
	(*par).Dpos[2] = new int [15];
	(*par).dimDposx = 3;
	(*par).dimDposy = (*par).numberGenes;
	
	(*par).Dpos[0][0] = 0;		
	(*par).Dpos[0][1] = -1;
	(*par).Dpos[0][2] = -1;
	(*par).Dpos[0][3] = 3;
	(*par).Dpos[0][4] = 4;
	(*par).Dpos[0][5] = 6;
	(*par).Dpos[0][6] = 7;
	(*par).Dpos[0][7] = 9;
	(*par).Dpos[0][8] = 10;
	(*par).Dpos[0][9] = 12;
	(*par).Dpos[0][10] = 14;
	(*par).Dpos[0][11] = 15;	
	(*par).Dpos[0][12] = 17;	
	(*par).Dpos[0][13] = -1;	
	(*par).Dpos[0][14] = -1;	

	(*par).Dpos[1][0] = 1;		
	(*par).Dpos[1][1] = -1;
	(*par).Dpos[1][2] = -1;
	(*par).Dpos[1][3] = -1;
	(*par).Dpos[1][4] = 5;
	(*par).Dpos[1][5] = -1;
	(*par).Dpos[1][6] = 8;
	(*par).Dpos[1][7] = -1;
	(*par).Dpos[1][8] = 11;
	(*par).Dpos[1][9] = 13;
	(*par).Dpos[1][10] = -1;
	(*par).Dpos[1][11] = 16;	
	(*par).Dpos[1][12] = -1;	
	(*par).Dpos[1][13] = -1;	
	(*par).Dpos[1][14] = -1;	
	
	(*par).Dpos[2][0] = 2;		
	(*par).Dpos[2][1] = -1;
	(*par).Dpos[2][2] = -1;
	(*par).Dpos[2][3] = -1;
	(*par).Dpos[2][4] = -1;
	(*par).Dpos[2][5] = -1;
	(*par).Dpos[2][6] = -1;
	(*par).Dpos[2][7] = -1;
	(*par).Dpos[2][8] = -1;
	(*par).Dpos[2][9] = -1;
	(*par).Dpos[2][10] = -1;
	(*par).Dpos[2][11] = -1;	
	(*par).Dpos[2][12] = -1;	
	(*par).Dpos[2][13] = -1;	
	(*par).Dpos[2][14] = -1;	
	
	
	
	// Setting up the Dneg index	
	(*par).Dneg = new int * [2];
	(*par).Dneg[0] = new int [15];
	(*par).Dneg[1] = new int [15];
	(*par).dimDnegx = 2;
	(*par).dimDnegy = (*par).numberGenes;
	
	(*par).Dneg[0][0] = -1;		
	(*par).Dneg[0][1] = 0;
	(*par).Dneg[0][2] = 2;
	(*par).Dneg[0][3] = -1;
	(*par).Dneg[0][4] = -1;
	(*par).Dneg[0][5] = 4;
	(*par).Dneg[0][6] = -1;
	(*par).Dneg[0][7] = 5;
	(*par).Dneg[0][8] = -1;
	(*par).Dneg[0][9] = -1;
	(*par).Dneg[0][10] = -1;
	(*par).Dneg[0][11] = -1;	
	(*par).Dneg[0][12] = -1;	
	(*par).Dneg[0][13] = 6;	
	(*par).Dneg[0][14] = 7;	

	(*par).Dneg[1][0] = -1;		
	(*par).Dneg[1][1] = 1;
	(*par).Dneg[1][2] = 3;
	(*par).Dneg[1][3] = -1;
	(*par).Dneg[1][4] = -1;
	(*par).Dneg[1][5] = -1;
	(*par).Dneg[1][6] = -1;
	(*par).Dneg[1][7] = -1;
	(*par).Dneg[1][8] = -1;
	(*par).Dneg[1][9] = -1;
	(*par).Dneg[1][10] = -1;
	(*par).Dneg[1][11] = -1;	
	(*par).Dneg[1][12] = -1;	
	(*par).Dneg[1][13] = -1;	
	(*par).Dneg[1][14] = -1;


	// Setting up the Upos index	
	(*par).Upos = new int * [3];
	(*par).Upos[0] = new int [15];
	(*par).Upos[1] = new int [15];
	(*par).Upos[2] = new int [15];
	(*par).dimUposx = 3;
	(*par).dimUposy = (*par).numberGenes;
	
	(*par).Upos[0][0] = 17;		
	(*par).Upos[0][1] = -1;
	(*par).Upos[0][2] = -1;
	(*par).Upos[0][3] = 0;
	(*par).Upos[0][4] = 1;
	(*par).Upos[0][5] = 4;
	(*par).Upos[0][6] = 2;
	(*par).Upos[0][7] = 14;
	(*par).Upos[0][8] = 9;
	(*par).Upos[0][9] = 5;
	(*par).Upos[0][10] = 11;
	(*par).Upos[0][11] = 8;	
	(*par).Upos[0][12] = -1;	
	(*par).Upos[0][13] = -1;	
	(*par).Upos[0][14] = -1;	

	(*par).Upos[1][0] = -1;		
	(*par).Upos[1][1] = -1;
	(*par).Upos[1][2] = -1;
	(*par).Upos[1][3] = -1;
	(*par).Upos[1][4] = 6;
	(*par).Upos[1][5] = 10;
	(*par).Upos[1][6] = 3;
	(*par).Upos[1][7] = 15;
	(*par).Upos[1][8] = 16;
	(*par).Upos[1][9] = -1;
	(*par).Upos[1][10] = 12;
	(*par).Upos[1][11] = 13;	
	(*par).Upos[1][12] = -1;	
	(*par).Upos[1][13] = -1;	
	(*par).Upos[1][14] = -1;	
	
	(*par).Upos[2][0] = -1;		
	(*par).Upos[2][1] = -1;
	(*par).Upos[2][2] = -1;
	(*par).Upos[2][3] = -1;
	(*par).Upos[2][4] = 7;
	(*par).Upos[2][5] = -1;
	(*par).Upos[2][6] = -1;
	(*par).Upos[2][7] = -1;
	(*par).Upos[2][8] = -1;
	(*par).Upos[2][9] = -1;
	(*par).Upos[2][10] = -1;
	(*par).Upos[2][11] = -1;	
	(*par).Upos[2][12] = -1;	
	(*par).Upos[2][13] = -1;	
	(*par).Upos[2][14] = -1;			


	// Setting up the Uneg index	
	(*par).Uneg = new int * [2];
	(*par).Uneg[0] = new int [15];
	(*par).Uneg[1] = new int [15];
	(*par).dimUnegx = 2;
	(*par).dimUnegy = (*par).numberGenes;
	
	(*par).Uneg[0][0] = -1;		
	(*par).Uneg[0][1] = 2;
	(*par).Uneg[0][2] = 7;
	(*par).Uneg[0][3] = -1;
	(*par).Uneg[0][4] = 0;
	(*par).Uneg[0][5] = 1;
	(*par).Uneg[0][6] = -1;
	(*par).Uneg[0][7] = 4;
	(*par).Uneg[0][8] = 3;
	(*par).Uneg[0][9] = -1;
	(*par).Uneg[0][10] = -1;
	(*par).Uneg[0][11] = -1;	
	(*par).Uneg[0][12] = -1;	
	(*par).Uneg[0][13] = -1;	
	(*par).Uneg[0][14] = -1;	

	(*par).Uneg[1][0] = -1;		
	(*par).Uneg[1][1] = 6;
	(*par).Uneg[1][2] = -1;
	(*par).Uneg[1][3] = -1;
	(*par).Uneg[1][4] = 5;
	(*par).Uneg[1][5] = -1;
	(*par).Uneg[1][6] = -1;
	(*par).Uneg[1][7] = -1;
	(*par).Uneg[1][8] = -1;
	(*par).Uneg[1][9] = -1;
	(*par).Uneg[1][10] = -1;
	(*par).Uneg[1][11] = -1;	
	(*par).Uneg[1][12] = -1;	
	(*par).Uneg[1][13] = -1;	
	(*par).Uneg[1][14] = -1;
	
	
	(*par).nDpos = new int [15];
	(*par).nDneg = new int [15];
	(*par).nUpos = new int [15];
	(*par).nUneg = new int [15];
	
	for (int i=0;i<(*par).numberGenes;i++)
	{
		(*par).nDpos[i] = 0;
		(*par).nDneg[i] = 0;
		(*par).nUpos[i] = 0;
		(*par).nUneg[i] = 0;
				
		for (int j=0;j<3;j++)
		{
			if ((*par).Upos[j][i]>=0)
				(*par).nUpos[i] ++;

			if ((*par).Dpos[j][i]>=0)
				(*par).nDpos[i] ++;
		}
		
		
		
		for (int j=0;j<2;j++)
		{

			if ((*par).Uneg[j][i]>=0)
				(*par).nUneg[i] ++;
				
			if ((*par).Dneg[j][i]>=0)
				(*par).nDneg[i] ++;

		}
		
		
		
		
	}
	
	
	// seting up the transcription rates
	
	// in this inital trial, we set the basal rate to be uniform:
	// 0.0 for genes which are only subject to activators (type 0)
	// 1.0 for genes which are only subject to repressors (type 1)
	// 0.5 for genes which are subject to both (type 2)
	// As for the diffTranscriptionRate we set it to be unirform as well:
	// 1.0 for activators acting on type 0 gene
	// -1.0 for repressors acting on type 1 gene
	// 0.5 for activators acting on type 2 gene
	// -0.5 for repressors acting on type 2 gene	
	
	double uniformRate = (*par).uniformTranscriptionRate;
	double basalRate = (*par).uniformBasal;
	double leakRate = (*par).uniformLeak;
	
	for (int i=0;i<(*par).numberGenes;i++)
	{
			
		if (((*par).nUpos[i]>0)&&((*par).nUneg[i]>0))
		{ 
			// this gene is regulated both negatively and positively
			(*par).basalTranscriptionRate[i] = basalRate;//0.05*uniformRate/10;
			(*par).geneType[i] = 2;
		}
		else if (((*par).nUpos[i]>0)&&((*par).nUneg[i]==0))
		{
			// this gene is only regulated positively
			(*par).basalTranscriptionRate[i] = leakRate;
			(*par).geneType[i] = 0;
		}
		else
		{
			// this gene is only regulated negatively
			(*par).basalTranscriptionRate[i] = uniformRate;
			(*par).geneType[i] = 1;
		}
			
	}
	
	
	// now we assign the diff transcription rates
	for (int i = 0; i<(*par).dimEposy;i++)
	{
		int indexTF = (*par).Epos[0][i];
		int indexGene = (*par).Epos[1][i];
		
		if ((*par).nUneg[indexGene]>0)
		{
			// indexGene are both regulated positively and negatively
			(*par).diffTranscriptionRate[indexTF][indexGene] = uniformRate-basalRate;
		}else{
			// indexGene is only regulated positively
			(*par).diffTranscriptionRate[indexTF][indexGene] = uniformRate-leakRate;
		}
	}
	
	for (int i = 0; i<(*par).dimEnegy;i++)
	{
		int indexTF = (*par).Eneg[0][i];
		int indexGene = (*par).Eneg[1][i];
		
		if ((*par).nUpos[indexGene]>0)
		{
			// indexGene are both regulated positively and negatively
			(*par).diffTranscriptionRate[indexTF][indexGene] = leakRate-basalRate;
		}else{
			// indexGene is only regulated positively
			(*par).diffTranscriptionRate[indexTF][indexGene] = leakRate-uniformRate;
		}
	}
	
	
	// finally, we set up koff and kon, assuming a uniform koff / kon distribution of all possible reactions
	double uniform_koff = (*par).uniformKoff;
	double uniform_kon =  (*par).uniformKon;

	for (int i = 0; i<(*par).dimEposy;i++)
	{
		int indexTF = (*par).Epos[0][i];
		int indexGene = (*par).Epos[1][i];
		
		(*par).koff[indexTF][indexGene] = uniform_koff;
		(*par).kon[indexTF][indexGene] = uniform_kon;
		
	}
	

	for (int i = 0; i<(*par).dimEnegy;i++)
	{
		int indexTF = (*par).Eneg[0][i];
		int indexGene = (*par).Eneg[1][i];
		
		(*par).koff[indexTF][indexGene] = uniform_koff;
		(*par).kon[indexTF][indexGene] = uniform_kon;
		
	}	  
	
	
	// degradation rate here
	
	for (int i=0;i<(*par).numberGenes;i++)
	{	
		(*par).gamma[i] = (*par).uniformGamma;	
	}
	
	
	// remove the dynamics of LIF, CH, PD
	(*par).basalTranscriptionRate[12] = (*par).gamma[12] * (*par).LIF;
	(*par).basalTranscriptionRate[13] = (*par).gamma[13] * (*par).CH;
	(*par).basalTranscriptionRate[14] = (*par).gamma[14] * (*par).PD;
	
	
	
}
