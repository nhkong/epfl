#!/usr/bin/python

import sys

firstPatientColumn = 44
numPatients = 176
threshold = 20

sys.stdout.write(sys.stdin.readline())

for line in sys.stdin:
        mutationSite = line.split("\t")
        patientIndex = firstPatientColumn

        #Write sites above threshold
        for patientIndex in range(firstPatientColumn, firstPatientColumn + numPatients):
            if (mutationSite[patientIndex] != "0/0" and int(mutationSite[patientIndex + numPatients]) > threshold):
                sys.stdout.write(line)
                break
        """
        #Write sites below threshold
        falsePositive = True
        for patientIndex in range(firstPatientColumn, firstPatientColumn + numPatients):
            if (mutationSite[patientIndex] != "0/0" and int(mutationSite[patientIndex + numPatients]) >= threshold):
                falsePositive = False
                break

        if (falsePositive):
            sys.stdout.write(line)
        """
