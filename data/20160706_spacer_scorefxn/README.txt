These simulations were run with a score function that would apply five 
different spacers (none, gfp, rfp, aavs, vegfa), get scores with each of them, 
then return the average score.

The purpose of this simulation was to see if there were positions that could 
make the switch more robust and less affected by the spacer sequence.

However, the resulting sequence logo doesn't look much different than the one 
for mh/7 (same sequence length) with the spacer ignored.  Since this simulation 
was also 5x slower, I don't think it's worth doing in the future.
