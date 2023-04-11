import sys
from pathlib import Path

import pandas as pd
from os import listdir
from os.path import isfile, join


def getBinName(f):
    return Path(f).stem

def main(args):
    print("args" + str(args))
    path =  "/work_beegfs/sukem127/clean/16ss/scafstats-statsfiles-unpaired-all/" #args[0]
    # Read files into a list of pandas  frames
    onlyfiles = list(map(lambda x: path + "/" + x, [f for f in listdir(path) if isfile(join(path, f))]))
    print(onlyfiles)
    dfs = list()
    # open files and read them into a df
    for f in onlyfiles:
        df = pd.read_csv(f, sep='\t')
        df = df[["#name", "unambiguousReads"]]
        df["#name"] = df["#name"].str.split(pat=" ", expand=True)[0]
        binName = getBinName(f) + ".fa"
        df = df.rename(columns={"#name" : "#name", "unambiguousReads" : binName})

        # Change name on the second colum to the binname so we can identify to which bin the amount of reads belongs
        # add df to df list(wasteful but whatever)
        dfs.append(df)

    df = dfs[0]
    for i in range(1, len(dfs) ):
        print("merging "   + str(i))
        df1 = dfs[i]

        df = df.merge(df1, on="#name")

    #print(dfs)
    df.to_csv("allstats-unambiguousReads.csv", sep = "\t", index=False)


if __name__ == '__main__':
    print("WD")
    args = sys.argv[1:]

    main(args)
    exit(0)
