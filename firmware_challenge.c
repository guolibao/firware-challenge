#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>


#define CYCLE 3
#define DATA_LENGTH 20

// Va = [156.63, 246.59, 294.72, 305.51, 300.66, 268.03, 204.18, 125.41, 42.954, -48.322, -154.08, -243.95, -293.12, -303.09, -297.98, -264.13, -202.1, -122.25, -39.893, 51.818];
float Va[] = {
    156.63, 246.59, 294.72, 305.51, 300.66,
    268.03, 204.18, 125.41, 42.954, -48.322,
    -154.08, -243.95, -293.12, -303.09, -297.98,
    -264.13, -202.1, -122.25, -39.893, 51.818
};

// Vb = [-308.4, -280.19, -240.66, -186.6, -99.744, -0.54547, 92.853, 181.46, 262.05, 312.39, 311.44, 283.76, 245.04, 188.62, 102.16, 2.9662, -89.395, -176.17, -259.16, -309.96];
float Vb[] = {
    -308.4, -280.19, -240.66, -186.6, -99.744,
    -0.54547, 92.853, 181.46, 262.05, 312.39,
    311.44, 283.76, 245.04, 188.62, 102.16,
    2.9662, -89.395, -176.17, -259.16, -309.96
};

// Vc = [156.11, 82.694, -21.783, -128.37, -213.06, -269.49, -309.58, -313.4, -273.73, -214.81, -154.29, -79.64, 24.679, 132.16, 216.63, 274.14, 311.11, 315.76, 276.27, 216.22];
float Vc[] = {
    156.11, 82.694, -21.783, -128.37, -213.06,
    -269.49, -309.58, -313.4, -273.73, -214.81,
    -154.29, -79.64, 24.679, 132.16, 216.63,
    274.14, 311.11, 315.76, 276.27, 216.22
};


typedef struct _DDATA {
    float* in_a;
    float* in_b;
    float* in_c;
    float F_est;
    float* Theta_est;
    float* Harmonics;
    float Ts;
    float Kc1; // Kc are controller gains
    float Kc2; // choose your controller and
    float Kc3; // gains accordingly to get satisfied result
}DDATA;

DDATA ddata;

typedef struct singleData {
    float data[2];
    int index[2];
    float interpIndex;
}singleData_t;


//! calculate the phase ange base on 0 crossing for Vb from all different sample points.
float phaseAngle_degree[DATA_LENGTH] = { 0 };


//! @ this function will do linear interpolation for the index for x 
//! linear interplotion equation y = y1 + (x-x1) * ((y2 - y1)/(x2-x1))
//! so if gain = (y2 -y1)/(x2 - x1)
//! then x = (y - y1)/ gain + x1
//! param[in] y2 - float value as above equation, y1 - value as above equation, y - input as above equation
//! @return float interpolated index
float linearInterpolateForIndex(float x1, float y1, float x2, float y2, float y)
{
    float gain = (y2 - y1) / (x2 - x1);
    float x = (y - y1) / gain + x1;
    return x;
}
// Function estimateFrequencyAndTheta calculates estimates
// for Theta_est and F_est on each data instance.
// Initially, the estimates may not be very accurate,
// but they shall improve with more input data.
// Thus, the size of the Estimated Theta and Freq arrays
// should match the input data size (DATALENGTH * CYCLE)
void estimateFrequencyAndTheta(DDATA* d, int dataSize) {
    // Implementation for estimating frequency and theta
        /* . . . */
    float a[DATA_LENGTH * CYCLE], b[DATA_LENGTH * CYCLE], c[DATA_LENGTH * CYCLE];
    for (int i=0; i < CYCLE; i++)
    { 
        for (int j = 0; j < DATA_LENGTH; j++)
        {
            int index = i * DATA_LENGTH + j;
            a[index] = d->in_a[j];
            b[index] = d->in_b[j];
            c[index] = d->in_c[j];
        }
    }
    singleData_t firstPair = { 0 };
    singleData_t secondPair = { 0 };
    bool firsPairFound = false;
    bool secondPairFound = false;
    float diff;
    // look for the 1st and 2nd rising edge to get the electrice period
    // we use pairs here, we assume, that no voltage data exactly lies on the zero cross
    // and one data is on the negative side, another is on the positive side. so we need to calculate the zero cross base on the two pair values 
    // by linear interpolation
    for (int i = 0; i < DATA_LENGTH * CYCLE; i++)
    {
        if (!firsPairFound)
        {
            if (b[i] < 0) // the rising edege, one point is at the negative side, another is at the postive side
            {
                firstPair.data[0] = b[i];
                firstPair.index[0] = i;
                if ((i + 1) < DATA_LENGTH * CYCLE)
                {
                    if (b[i + 1] > 0)
                    {
                        firstPair.data[1] = b[i + 1];
                        firstPair.index[1] = i + 1;
                        firsPairFound = true;
                    }
                }
            }
        }
        else
        {
           // first pair is found, look for the 
            if (!secondPairFound)
            {
                if (b[i] < 0) // the rising edege, one point is at the negative side, another is at the postive side
                {
                    secondPair.data[0] = b[i];
                    secondPair.index[0] = i;
                    if ((i + 1) < DATA_LENGTH * CYCLE)
                    {
                        if (b[i + 1] > 0)
                        {
                            secondPair.data[1] = b[i + 1];
                            secondPair.index[1] = i + 1;
                            secondPairFound = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    if (firsPairFound && secondPairFound) // after found, then calculate the two consecutive zero cross index
    {
        // do the interplotion to get the 0 cross index, assume it is linear
        // x1 = tempNegtive[0].index, y1 = tempNegtive[0].data 
        // x2 = tempPositive[0].index, y2 = tempPositive[0].data
        // x = ?                      , y = 0
        // linear interplotion equation y = y1 + (x-x1) * ((y2 - y1)/(x2-x1))
        // so if gain = (y2 -y1)/(x2 - x1)
        // then x = (y - y1)/ gain + x1
        printf("first pair: x1 = %d, y1=%.2f, x2=%d, y2=%.2f\r\n", firstPair.index[0], firstPair.data[0], firstPair.index[1], firstPair.data[1]);
        printf("2nd pair: x1 = %d, y1=%.2f, x2=%d, y2=%.2f\r\n", secondPair.index[0], secondPair.data[0], secondPair.index[1], secondPair.data[1]);

        firstPair.interpIndex = linearInterpolateForIndex(firstPair.index[0], firstPair.data[0], firstPair.index[1], firstPair.data[1], 0.0f);
        printf("first interpolated zero cross index = %f\r\n", firstPair.interpIndex);

        secondPair.interpIndex = linearInterpolateForIndex(secondPair.index[0], secondPair.data[0], secondPair.index[1], secondPair.data[1], 0.0f);
        printf("second interpolated zero cross index = %f\r\n", secondPair.interpIndex);
    }
    // after found the zero cross index, we need to calculate electric period in seconds
    float electricPeriod_s = (secondPair.interpIndex - firstPair.interpIndex) * d->Ts;
    printf("one electric period = %f seconds\r\n", electricPeriod_s);
    // convert the period to Hz
    d->F_est = 1.0f / electricPeriod_s;
    printf("Frequency = %.2f Hz\r\n", d->F_est);

    // now we know, one electric period is electricPeriod_s
    // e.g. 360(or 2\pi), it will take such amount of time (electricPeriod_s)
    // so if it rotates 1 degree, it takes (electricPeriod_s/360) time in s
    // one Ts period, it will travel (electricPeriod_s/360) * Ts
    float oneTsRotation_degree = (360.0f/electricPeriod_s) * d->Ts;

    // we take Vb as an example to calculate all the phase angles for different sample points
    // because the 0 crossing calculation is base on the Vb sample data.
    for (int i = 0; i < DATA_LENGTH; i++)
    {
        if (i < firstPair.interpIndex)
        {
            // the nearer to the interpIndex, the nearer to 360 or 0 degrees
            phaseAngle_degree[i] = 360 - (firstPair.interpIndex - i) * oneTsRotation_degree;
        }
        else if (i == firstPair.interpIndex)
        {
            phaseAngle_degree[i] = 0;
        }
        else
        {
            // i > than interpIndex, the angle maybe greater than 360, but we make the angle within 360 degrees by using modular
            float degrees = (i - firstPair.interpIndex) * oneTsRotation_degree;
            phaseAngle_degree[i] = (int)degrees % 360;
        }
    }
    d->Theta_est = phaseAngle_degree;
    for (int i = 0; i < DATA_LENGTH; i++)
    {
        printf("phaseAngle_degree[%d]=%7.2f\r\n", i, d->Theta_est[i]);
    }
}


// Function getHarmonicAmplitudes calculates the amplitudes of the 1st to 5th
// The output should consist of 5 data points, each representing the amplitud
// Additionally, you can use these 5 data points to calculate the Total Harmo
void getHarmonicAmplitudes(DDATA* d, int dataSize) {
    // Implementation for getting harmonic amplitudes
        /* . . . */
    // use fft to get amplitude for each homonics
}

//! @brief Assume the function Va, Vb, or Vc data points is taken at time 0s
//! This function is to find the first rising edge,
//! two situations, 
//! (1). two points, the first point is negative and the next to this point is positive.
//!      we need to get these two point and linear intepolate to get the 0 crossing time
//! (2). one point, the y value is exactly 0.
//! @return the time in seconds when the first zero cross happens

int main()
{
    int i = 0;

    ddata.in_a = Va;
    ddata.in_b = Vb;
    ddata.in_c = Vc;
    ddata.Ts = 0.001;

    //for (i = 0; i < DATA_LENGTH * CYCLE; i++)
    {
        estimateFrequencyAndTheta(&ddata, DATA_LENGTH * CYCLE);
    }
    getHarmonicAmplitudes(&ddata, DATA_LENGTH * CYCLE);
    return 0;
}


