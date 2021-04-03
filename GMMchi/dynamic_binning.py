import numpy as np
import pandas as pd



def dynamic_binning(observed, binedges, threshold=5, filt=0, final=True):
    continued = 1
    x = 0  # dynamically update range
    orgi = observed.copy()

    while x < len(observed):
        try:
            restart = True  # loop until larger than threshold
            while restart:

                if observed[x] < threshold:
                    observed[x] = observed[x]+observed[x+1]
                    observed = np.delete(observed, x+1)
                    binedges = np.delete(binedges, x+1)
                else:
                    restart = False
            x += 1
        except:  # sometimes you will get bins with >threshold in the very last bin
            if observed[-1] < threshold:

                observed[-1] = observed[-1]+observed[-2]
                observed = np.delete(observed, -2)
                binedges = np.delete(binedges, -2)

    if len(orgi) == len(observed):
        continued = 0

    if final == False:
        largedata = np.arange(len(observed))[np.in1d(
            observed, orgi)]  # get overlapping index

        try:
            # check whether it is spanning the entire dataset
            if min(largedata) != 0 or max(largedata)+1 != len(observed):

                # expand the space where there is no tail problem to the number of bins using mann and wald and leave the rest
                if len(largedata) == 1:
                    # mann and wald
                    nums = int(1.88*observed[largedata[0]]**(2/5))
                    x = largedata[0]
                    newbins = np.linspace(
                        binedges[largedata[0]], binedges[largedata[0]+1], num=nums)
                    binedges = np.delete(
                        binedges, (largedata[0], largedata[0]+2))
                    binedges = np.sort(np.concatenate((binedges, newbins)))
                else:
                    # mann and wald
                    nums = int(
                        1.88*sum(observed[min(largedata):max(largedata)+1])**(2/5))
                    newbins = np.linspace(
                        binedges[min(largedata)], binedges[max(largedata)+1], num=nums)
                    binedges = np.delete(binedges, np.arange(
                        min(largedata), max(largedata)+2))
                    binedges = np.sort(np.concatenate((binedges, newbins)))
        except:
            pass  # if they are the same during the middle stage

    return observed, binedges, continued
