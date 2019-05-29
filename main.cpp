

/*---------------------------------------------------------------------------------*/

/*Dragged Overdamped Particle in a Hear Reservior*/

/*---------------------------------------------------------------------------------*/


#include  <fstream>
#include  <iostream>
#include <string>
#include  <cmath>
#include  <random>
#include <cfenv>
#include <climits>
#include <iomanip>

#include "nr3.h"
#include "gamma.h"
#include "gaussj.h"
#include "fitmrq.h"
#include "fit_examples.h"

#define PI     3.1415925358793
#define KB     1.38064852E-23
#define NN     1600000
#define NS     1600000
#define NRANMX 3200000
#define BIN    1000
#define NVSTEP 1
#define KBT    300

static double RAN[NRANMX];
static double position[NN];
static double EquPosition[NN];

void RANDOMCAL();

int main(void)
{
    //initializing the parameters
    double VEL[NVSTEP] = {};
    double deltaV = 10.0;
    double k = 7.0, eta = 1E-1, mass = 1.0;
    int ocillationNumber =5000;
    int simNumber = 100;  //This number as in (simTime / simNumber) define how long we should wait till the systems goes to equilibrium
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int NumberOfSteps = 20000;          // Number of iteration per oscillation
    int NumofStepsPerSample = 500;      //data are sampled once every 500 time steps
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    double period = (2.0 * PI) * sqrt(mass / k);
    double timeStep = period / NumberOfSteps;
    double simTime = ocillationNumber * period;
    
    double maxPosition = -100.0, minPosition = 100.0;
    

    
    
    
    std::ofstream  outputFile, outputFile2, outputFile3, outputFile4;
    // outputFile.open ("ParticlePosition.dat");
    // outputFile2.open("EquilibriumPosition.dat");
    // outputFile3.open("Probabilityfunction.dat");
    outputFile4.open("Force_Temperture.txt");
    
    std::cout << "What kind of system do you want to consider? "<<"\n"<<"input `1' for damped  and `2' for overdamped"<< "\n";
    //std::cin >> method;
    int method = 1;
    
    for (int m = 0; m < NVSTEP; m++)
    {
        
        std::cout<< "m = "<< m << "\n";
        std::cout<<"VEL = "<<VEL[m]<<"\n";
        
        outputFile.open("ParticlePosition.dat");
        outputFile2.open("EquilibriumPosition.dat");
        outputFile3.open("Probabilityfunction.dat");
        
        double ElapseTime = 0.0;
        
        int j = 0, l = 0;
        
        VEL[m]= m*deltaV;

        int NRAN = 0;
        RANDOMCAL();
        
        double v = 0.0;
        double x = 200;
        int ntime = 0;
        switch (method)
        {
            case 1:  //This is the program for the damped oscillator
            {
                //setting up the equilibrium
                while (ElapseTime <= simTime / double(simNumber))
                {
                    ElapseTime += timeStep;
                    ntime++;
                    v = v - timeStep*(k / mass)*(x)-timeStep*(eta / mass)*v + timeStep*(eta / mass)*VEL[m] + sqrt(4.0 * double(KBT)/eta*timeStep)*RAN[NRAN];
                    x = x + timeStep*v;
                    NRAN++;
                    
                    if (NRAN >= NRANMX / 2)
                    {
                        RANDOMCAL();
                        NRAN = 0;
                    }
                    
                    if ((ntime % NumofStepsPerSample) == 0)
                    {
                        
                        position[j] = x;
                        
                        //std::cout << position[j] << "\n";
                        //std::cin.get();
                        //ElapseTime = ntime*timeStep;
                        outputFile << std::left << std::fixed << std::setw(32) << ElapseTime << std::setfill(' ') << std::setw(16) << std::setprecision(16) << position[j] << "\n";
                        j++;
                    }
                }
                
                //the main loop
                ntime = 1;
                while (ElapseTime > simTime / double(simNumber) && ElapseTime <= simTime)
                {
                    ntime++;
                    ElapseTime += timeStep;
                    v = v - timeStep*(k / mass)*(x)-timeStep*(eta / mass)*v + timeStep*(eta / mass)*VEL[m] + sqrt(4.0 * double(KBT)/(eta*timeStep))*RAN[NRAN];
                    x = x + timeStep*v;
                    NRAN++;
                    
                    if (NRAN >= NRANMX / 2)
                    {
                        RANDOMCAL();
                        NRAN = 1;
                    }
                    
                    
                    if ((ntime % NumofStepsPerSample) == 0)
                    {
                        //std::cout << ElapseTime;
                        EquPosition[l] = x;
                       // outputFile2 << std::left << std::fixed << std::setw(32) << ElapseTime << std::setfill(' ') << std::setw(16) << std::setprecision(16) << EquPosition[l] << "\n";
                        if (maxPosition < EquPosition[l])  maxPosition = EquPosition[l];
                        if (minPosition > EquPosition[l])  minPosition = EquPosition[l];
                        l++;
                    }
                }
                
                //Histogram of the results
                double np[BIN];
                double bin_position[BIN];
                double delta = (maxPosition - minPosition) / BIN;
                for (int i = 0; i < BIN; i++)
                {
                    np[i] = 0.0;
                    bin_position[i] = 0.0;
                }
                for (int i = 0; i < BIN; i++)
                {
                bin_position[i] = minPosition + (2.0 * double(i) + 1) / 2.0 * delta;
                
                    for (int j = 0; j < l; j++)
                    {
                        if (EquPosition[j] <= minPosition + i*delta && EquPosition[j] > minPosition + (i - 1)*delta)
                        {
                            np[i] = np[i] + 1;
                        }
                    }
                }

                
                for (int i = 0; i < BIN; i++)
                {
                    outputFile3 << std::left << std::fixed << std::setw(32) << minPosition + (2.0 * double(i) + 1) / 2.0 * delta << ' ' << std::setw(16) << np[i] << "\n";
                }
                //calculating the guess parameters in the gaussian fit function.
                double S[BIN], T[BIN];
                double M[2][2];
                S[0] = 0.0;
                T[0] = 0.0;
                for (int k = 1; k < BIN; k++)
                {
                    S[k] = S[k - 1] + .5*(np[k] + np[k - 1])*(bin_position[k] - bin_position[k - 1]);
                    T[k] = T[k - 1] + .5*(np[k] * bin_position[k] + np[k - 1] * bin_position[k - 1])*(bin_position[k] - bin_position[k - 1]);
                }
                double m1 = 0.0, m2 = 0.0,  m4 = 0.0;
                double g1 = 0.0, g2 = 0.0;
                
                for (int k = 0; k < BIN; k++)
                {
                    m1 += S[k] * S[k];
                    m2 += S[k] * T[k];
                    m4 += T[k] * T[k];
                    
                    g1 += (np[k] - np[0])*S[k];
                    g2 += (np[k] - np[0])*T[k];
                    
                }
                //finding the inverse matrix
                double delta2 = (m1*m4 - m2*m2);
                M[0][0] =  m4 / delta2;
                M[0][1] = -m2 / delta2;
                M[1][0] = -m2 / delta2;
                M[1][1] =  m1 / delta2;
                
                //Finding Force and Temperture
                
                double a = M[0][0] * g1 + M[0][1] * g2;
                double b = M[0][1] * g1 + M[1][1] * g2;
                
                double force = (-a / b)*k;
                double temperture = k/2.0*(-2.0 / b);
                
                //Finding the amplitude
                double sum1=0, sum2=0;
                for(int k=0; k<BIN;k++ )
                {
                    sum1 += np[k];
                    sum2 += exp(-pow(bin_position[k]-(-a/b), 2.0)/(-2.0 / b));
                }
                double amplitude = sum1/sum2;
                
                std::cout   <<"VEL[m]="<<VEL[m]<<"    "<< "force="<< force << "   " <<"Temperture="<<temperture <<"     "<<"amplitude="<< amplitude<<"\n";
                outputFile4 <<VEL[m]<<"    "<< force << "   " << "      "<<temperture <<"     "<< amplitude<<"\n";
                
                /**************************************************************************/
                //fitting the Data to a Gaussian equation using Levenberg-Marquardt Method//
                
                VecDoub      xx(BIN);
                VecDoub      y(BIN);
                VecDoub      sd(BIN);
                
                const Doub guess_d[]  = {amplitude, force/k, sqrt(temperture*(2.0/k))};
                VecDoub  guess(3, guess_d);
                
                for(int k=0; k<BIN; k++)
                {
                    xx[k]  =  bin_position[k];
                    y[k]   =  np[k];
                    sd[k]  =  0.001;
                }
                
                Fitmrq mymrq(xx, y, sd, guess, fgauss); // Use default tolerance 1.0e-3
                mymrq.fit();
                
                //Writing down the fitting results
                cout << "Initial guess: ";
                for (int i = 0; i < guess.size(); i++) {
                    cout << guess[i] << "  ";
                }
                cout << endl;
                cout << scientific << setprecision(5);
                cout << setw(10) << "Guessed"
                << setw(18) << "Fitted"
                << "      "
                << setw(13) << "Err. Est."
                << endl;
                for (int i = 0; i < guess.size(); i++)
                {
                    cout << setw(13) << guess[i]
                    << setw(16) << mymrq.a[i]
                    << "  +-  " << sqrt(mymrq.covar[i][i])
                    << endl;
                }
                cout << endl;
                cout << fixed;
                cout << "Number of points = " << BIN
                << ", Chi-squared = " << mymrq.chisq << endl;

            }		
        }
        outputFile.close();
        outputFile2.close();
        outputFile3.close();
        
    }
    return 0;
}


void RANDOMCAL()
{
    int i;
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    
    for (i = 1; i <= NRANMX; i++)
    {
        std::normal_distribution<double>   distribution(0, 0.71);
        RAN[i] = distribution(generator) ;
    }
    
}

