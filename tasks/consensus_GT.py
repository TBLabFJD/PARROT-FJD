
# Author: Gonzalo Núñez Moreno
# Date: 07/03/2023
# Script to find the consensus Genotype between callers
import pandas as pd
import sys

df = pd.read_table(sys.argv[1], sep="\t", header=None)
# df = pd.read_table("prueba2_GT.txt", sep="\t", header=None)

# Each caller has a differnet weigth so that in case of tie, the order of priorization are: dragen (2nd, deepvariant (1st), and GATK (3rd)
df_new = pd.concat([df[0]] * 3 + [df[1]] * 4 + [df[2]] * 2 , axis=1, ignore_index=True)

# print(df_new.mode(axis=1)[0])

# Write consensus Genotype
df_new.mode(axis=1).to_csv(sys.argv[2], sep='\t', header=False, index=False)


# Calculate and write discordances in genotype
discordance_list=list(map(lambda x: len(set(filter(lambda y: y == y , x))) - 1,df.values))
with open(sys.argv[3], 'w') as fp:
    for item in discordance_list:
        # write each item on a new line
        fp.write("%s\n" % item)



