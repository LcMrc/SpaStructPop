#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Initialization(int WildTypeVector[], int MutantVector[], int InitialNumberOfWildTypes, int InitialNumberOfMutants, double *Time)
{
	WildTypeVector[0] = InitialNumberOfWildTypes;
	WildTypeVector[1] = InitialNumberOfWildTypes;
	WildTypeVector[2] = InitialNumberOfWildTypes;
	WildTypeVector[3] = InitialNumberOfWildTypes;
	WildTypeVector[4] = 0;
	MutantVector[0] = 0;
	MutantVector[1] = 0;
	MutantVector[2] = 0;
	MutantVector[3] = 0;
	MutantVector[4] = InitialNumberOfMutants;	
	*Time = 0;
}

int ComputeSum(int PopulationVector[], int NumberOfDemes, int *Sum)
{

	int i, q = 0;

	for (i = 0; i < NumberOfDemes; i++)
	{ 
		q = q+PopulationVector[i];
	}

	*Sum = q;

}

double ComputeTotalDivisionRate(double Fitness, int WildTypeVector[], int MutantVector[], int CarryingCapacity, int NumberOfDemes)
{
	double Sum = 0;
	int i;

	for (i = 0; i < NumberOfDemes; i++)
	{
		Sum = Sum+Fitness*(1-(double)(WildTypeVector[i]+MutantVector[i])/CarryingCapacity)*WildTypeVector[i];
	}

	return Sum;
}

double ComputeTotalDeathRate(double DeathRate, int PopulationVector[], int NumberOfDemes)
{
	double Sum = 0;
	int i;

	for (i = 0; i < NumberOfDemes; i++)
	{
		Sum = Sum+DeathRate*PopulationVector[i];
	}

	return Sum;
}

double ComputeTotalMigrationRate(int WildTypeVector[], int MutantVector[], double FromCenterToLeafMigrationRate, double FromLeafToCenterMigrationRate)
{
	double Sum = 0;

	Sum = 4*FromCenterToLeafMigrationRate*(WildTypeVector[0]+MutantVector[0])+FromLeafToCenterMigrationRate*(WildTypeVector[1]+MutantVector[1]+WildTypeVector[2]+MutantVector[2]+WildTypeVector[3]+MutantVector[3]+WildTypeVector[4]+MutantVector[4]);

	return Sum;
}

void TotalTransitionRate(double WildTypeFitness, double WildTypeDeathRate, int WildTypeVector[], double MutantFitness, double MutantDeathRate, int MutantVector[], int CarryingCapacity, int NumberOfDemes, double FromLeafToCenterMigrationRate, double FromCenterToLeafMigrationRate, double *DivisionOrDeathRateVector)
{
	DivisionOrDeathRateVector[0] = ComputeTotalDivisionRate(WildTypeFitness, WildTypeVector, MutantVector, CarryingCapacity, NumberOfDemes);
	DivisionOrDeathRateVector[1] = ComputeTotalDeathRate(WildTypeDeathRate, WildTypeVector, NumberOfDemes);
	DivisionOrDeathRateVector[2] = ComputeTotalDivisionRate(MutantFitness, MutantVector, WildTypeVector, CarryingCapacity, NumberOfDemes);
	DivisionOrDeathRateVector[3] = ComputeTotalDeathRate(MutantDeathRate, MutantVector, NumberOfDemes);
	DivisionOrDeathRateVector[4] = ComputeTotalMigrationRate(WildTypeVector, MutantVector, FromCenterToLeafMigrationRate, FromLeafToCenterMigrationRate);
}

void ComputeCumulativeSum(double *DivisionOrDeathRateVector, int NumberOfReactions) 
{
    if(NumberOfReactions <= 0) return;
    ComputeCumulativeSum(DivisionOrDeathRateVector, NumberOfReactions-1);
    DivisionOrDeathRateVector[NumberOfReactions] += DivisionOrDeathRateVector[NumberOfReactions-1];
}

double GenerateRandomNumber()
{
    return (double)rand()/(double)RAND_MAX;
}

int SamplingReaction(double *DivisionOrDeathRateVector, int NumberOfReactions) 
{
	int i = 0;
	double RandomNumber1, RandomNumber2;
	
	RandomNumber1 = GenerateRandomNumber();
	RandomNumber2 = RandomNumber1*DivisionOrDeathRateVector[NumberOfReactions-1];

	while (DivisionOrDeathRateVector[i] < RandomNumber2)
	{
		i++;
	}

	return i;
}

double UpdateTime(double *Time, double TotalTransitionRate)
{
	double RandomNumber;
	RandomNumber = GenerateRandomNumber();
	*Time += 1/TotalTransitionRate*log(1/RandomNumber);
}

void ComputeDivisionOrDeathRate(double *DivisionOrDeathRateVector, int *WildTypeVector, int *MutantVector, int CarryingCapacity, int NumberOfDemes, int Reaction)
{
	int i;

	if (Reaction == 0)
	{
		for (i = 0; i < NumberOfDemes; i++)
		{
			DivisionOrDeathRateVector[i] = (1-(double)(WildTypeVector[i]+MutantVector[i])/CarryingCapacity)*WildTypeVector[i];
		}
	}
	else if (Reaction == 1)
	{
		for (i = 0; i < NumberOfDemes; i++)
		{
			DivisionOrDeathRateVector[i] = WildTypeVector[i];
		}
	}
	else if (Reaction == 2)
	{
		for (i = 0; i < NumberOfDemes; i++)
		{
			DivisionOrDeathRateVector[i] = (1-(double)(WildTypeVector[i]+MutantVector[i])/CarryingCapacity)*MutantVector[i];
		}
	}
	else
	{
		for (i = 0; i < NumberOfDemes; i++)
		{
			DivisionOrDeathRateVector[i] = MutantVector[i];
		}
	}
}

void PartialTransitionRateTer(double *DivisionOrDeathRateVector, int *WildTypeVector, int *MutantVector, double FromCenterToLeafMigrationRate, double FromLeafToCenterMigrationRate)
{
	DivisionOrDeathRateVector[0] = FromCenterToLeafMigrationRate*(WildTypeVector[0]+MutantVector[0]);
	DivisionOrDeathRateVector[1] = FromCenterToLeafMigrationRate*(WildTypeVector[0]+MutantVector[0]);
	DivisionOrDeathRateVector[2] = FromCenterToLeafMigrationRate*(WildTypeVector[0]+MutantVector[0]);
	DivisionOrDeathRateVector[3] = FromCenterToLeafMigrationRate*(WildTypeVector[0]+MutantVector[0]);
	DivisionOrDeathRateVector[4] = FromLeafToCenterMigrationRate*(WildTypeVector[1]+MutantVector[1]);
	DivisionOrDeathRateVector[5] = FromLeafToCenterMigrationRate*(WildTypeVector[2]+MutantVector[2]);
	DivisionOrDeathRateVector[6] = FromLeafToCenterMigrationRate*(WildTypeVector[3]+MutantVector[3]);
	DivisionOrDeathRateVector[7] = FromLeafToCenterMigrationRate*(WildTypeVector[4]+MutantVector[4]);	
}

void ExecuteDivisionOrDeath(int *WildTypeVector, int *MutantVector, int Reaction, int TargetedDeme)
{
	if (Reaction == 0)
	{
		WildTypeVector[TargetedDeme] += 1;
	}
	else if (Reaction == 1)
	{
		WildTypeVector[TargetedDeme] -= 1;
	}
	else if (Reaction == 2)
	{
		MutantVector[TargetedDeme] += 1;
	}
	else
	{
		MutantVector[TargetedDeme] -= 1;
	}
}

void ExecuteMigration(int *WildTypeVector, int *MutantVector, int TargetedDeme)
{
	int DepartureDeme, ArrivalDeme;
	double RandomNumber;

	if (TargetedDeme == 0)
	{
		DepartureDeme = 0;
		ArrivalDeme = 1;
	}
	else if (TargetedDeme == 1)
	{
		DepartureDeme = 0;
		ArrivalDeme = 2;
	}
	else if (TargetedDeme == 2)
	{
		DepartureDeme = 0;
		ArrivalDeme = 3;
	}
	else if (TargetedDeme == 3)
	{
		DepartureDeme = 0;
		ArrivalDeme = 4;
	}
	else if (TargetedDeme == 4)
	{
		DepartureDeme = 1;
		ArrivalDeme = 0;
	}
	else if (TargetedDeme == 5)
	{
		DepartureDeme = 2;
		ArrivalDeme = 0;
	}
	else if (TargetedDeme == 6)
	{
		DepartureDeme = 3;
		ArrivalDeme = 0;
	}
	else if (TargetedDeme == 7)
	{
		DepartureDeme = 4;
		ArrivalDeme = 0;
	}

	RandomNumber = GenerateRandomNumber();

	if (RandomNumber <= WildTypeVector[DepartureDeme]/(WildTypeVector[DepartureDeme]+MutantVector[DepartureDeme]))
	{
		WildTypeVector[DepartureDeme] -= 1;
		WildTypeVector[ArrivalDeme] += 1;
	}
	else
	{
		MutantVector[DepartureDeme] -= 1;
		MutantVector[ArrivalDeme] += 1;
	}
}

int main()
{

	const double WildTypeFitness = 1.0, WildTypeDeathRate = 0.1, MutantFitness = 1.01, MutantDeathRate = 0.1, FromLeafToCenterMigrationRate = 0.000001, FromCenterToLeafMigrationRate = 0.000001;
	const int CarryingCapacity = 100, InitialNumberOfWildTypes = 90, InitialNumberOfMutants = 90, NumberOfDemes = 5, NumberOfRealizations = 1000;
	double Time = 0, TransitionRateVector[5], DivisionOrDeathRateVector[5], MigrationRateVector[8], FinalTimeVector[1000];
	int WildTypeVector[5] = {0, InitialNumberOfWildTypes, InitialNumberOfWildTypes, InitialNumberOfWildTypes, InitialNumberOfWildTypes}, MutantVector[5] = {InitialNumberOfMutants, 0, 0, 0, 0}, k, i, j, Reaction, TargetedDeme, TotalNumberOfWildTypes = 0, TotalNumberOfMutants = 0, FinalMutantNumberVector[1000];

	for (k = 0; k < NumberOfRealizations; k++)
	{
		Initialization(WildTypeVector, MutantVector, InitialNumberOfWildTypes, InitialNumberOfMutants, &Time);
		ComputeSum(WildTypeVector, NumberOfDemes, &TotalNumberOfWildTypes);
		ComputeSum(MutantVector, NumberOfDemes, &TotalNumberOfMutants);

		while (TotalNumberOfWildTypes != 0 && TotalNumberOfMutants != 0)
		{
			TotalTransitionRate(WildTypeFitness, WildTypeDeathRate, WildTypeVector, MutantFitness, MutantDeathRate, MutantVector, CarryingCapacity, NumberOfDemes, FromLeafToCenterMigrationRate, FromCenterToLeafMigrationRate, TransitionRateVector);
			ComputeCumulativeSum(TransitionRateVector, NumberOfDemes-1);
			Reaction = SamplingReaction(TransitionRateVector, NumberOfDemes); 

			if (Reaction < 4)
			{
				ComputeDivisionOrDeathRate(DivisionOrDeathRateVector, WildTypeVector, MutantVector, CarryingCapacity, NumberOfDemes, Reaction);
				ComputeCumulativeSum(DivisionOrDeathRateVector, NumberOfDemes-1);
				TargetedDeme = SamplingReaction(DivisionOrDeathRateVector, NumberOfDemes);
				ExecuteDivisionOrDeath(WildTypeVector, MutantVector, Reaction, TargetedDeme);
			}
			else
			{
				PartialTransitionRateTer(MigrationRateVector, WildTypeVector, MutantVector, FromCenterToLeafMigrationRate, FromLeafToCenterMigrationRate);
				ComputeCumulativeSum(MigrationRateVector, 8-1);
				TargetedDeme = SamplingReaction(MigrationRateVector, 8);
				ExecuteMigration(WildTypeVector, MutantVector, TargetedDeme);
			}

			UpdateTime(&Time, TransitionRateVector[NumberOfDemes-1]);

			ComputeSum(WildTypeVector, NumberOfDemes, &TotalNumberOfWildTypes);
			ComputeSum(MutantVector, NumberOfDemes, &TotalNumberOfMutants);	 
		}
		
		FinalMutantNumberVector[k] = MutantVector[4];
		FinalTimeVector[k] = Time;
	}

	FILE *fp1;
	fp1 = fopen("FinalMutantNumberVector.txt", "w");
	for (i = 0; i < NumberOfRealizations; i++)
	{
		fprintf(fp1, "%d\n", FinalMutantNumberVector[i]);
	}
	fclose(fp1); 

	FILE *fp2;
	fp2 = fopen("FinalTimeVector.txt", "w");
	for (i = 0; i < NumberOfRealizations; i++)
	{
		fprintf(fp2, "%f\n", FinalTimeVector[i]);
	}
	fclose(fp2);

	return 0;
} 
