for (unsigned int i = 0; i < allLut.size(); i++)
    if (allLut[i].gateName.compare(gateName) == 0)
    {
        unsigned int j = 0;

        // slew is greater than last value in index1
        if (tau > allLut[i].inSlew[allLut[i].inSlew.size() - 1])
        {
            p = allLut[i].inSlew.size() - 2;
            q = allLut[i].inSlew.size() - 1;
        }

        // slew is smaller than first value of index1
        else if (tau < allLut[i].inSlew[0])
        {
            p = 0;
            q = 1;
        }

        // slew is with the range of values of index1
        else
        {
            while (allLut[i].inSlew[j++] <= tau);
                p = j - 2;
                q = j - 1;
        }

        j = 0;
        // cap is greater than last value in index2
        if (c > allLut[i].outLoad[allLut[i].outLoad.size() - 1])
        {
            r = allLut[i].outLoad.size() - 2;
            s = allLut[i].outLoad.size() - 1;
        }

        // cap is lesser than first value in index2
        else if (c < allLut[i].outLoad[0])
        {
            r = 0;
            s = 1;
        }

        // cap is with the range of index2
        else
        {
            while (allLut[i].outLoad[j++] <= c);
                r = j - 2;
                s = j - 1;
        }

        // all the values required for interpolation
        double tau1, tau2, c1, c2, v11, v12, v21, v22, v;
        tau1 = allLut[i].inSlew[p];
        tau2 = allLut[i].inSlew[q];
        c1 = allLut[i].outLoad[r];
        c2 = allLut[i].outLoad[s];
        v11 = allLut[i].delay[p][r];
        v12 = allLut[i].delay[p][s];
        v21 = allLut[i].delay[q][r];
        v22 = allLut[i].delay[q][s];

        // delay interpolation
        v = (v11*(c2 - c)*(tau2 - tau) + v12*(c - c1)*(tau2 - tau) + v21*(c2 - c)*(tau - tau1) + v22*(c - c1)*(tau - tau1)) / ((c2 - c1) * (tau2 - tau1));
        double vDelay = v;

        v11 = allLut[i].outSlew[p][r];
        v12 = allLut[i].outSlew[p][s];
        v21 = allLut[i].outSlew[q][r];
        v22 = allLut[i].outSlew[q][s];

        // slew interpolation
        v = (v11*(c2 - c)*(tau2 - tau) + v12*(c - c1)*(tau2 - tau) + v21*(c2 - c)*(tau - tau1) + v22*(c - c1)*(tau - tau1)) / ((c2 - c1) * (tau2 - tau1));
        double vSlew = v;

        // Make sure delay is non-negative and slew is greater than 2ps
        if (vSlew < 0.002) vSlew = 0.002;
        if (vDelay < 0) vDelay = 0;

        pair<double, double> valuePair;

        // Adjust delay and slew for n>2 input gates
        if (fanIn > 2)
            vDelay = vDelay * fanIn / 2;
        
        valuePair.first = vDelay;
        
        if (fanIn > 2)
            vSlew = vSlew * fanIn / 2;
        
        valuePair.second = vSlew;

        return(valuePair);
}

// Cannot interpolate because the gate does not exist in LUT
cout << "ERROR: Cannot find out delay for " << gateName << " in liberty file\n";
exit(0);
}