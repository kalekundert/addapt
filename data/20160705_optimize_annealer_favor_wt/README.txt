First, I chose a cycle length.  I ran simulations with cycle lengths of 
N=100,200,300,400,500.  For each simulation the high temperature was 5 and the 
low temperature was 0.  I chose the high temperature based on "auto-scaling" 
simulations I ran while tuning the score function.

I looked at the score trajectory for each and judged them based on how high the 
peaks were and how long the peaks lasted for.  For N=500 and N=600, the peaks 
seemed to last for a substantial fraction of the total simulation time, which 
was wasteful.

To corroborate this, I looked at the auto-correlation time graphs.  However, 
contrary to what I expected, all the correlation times were between 200-250 
steps except for N=200, which was maybe as low as 150.  The correlation time 
did not seem to increase as the schedule got slower.

I chose to move forward with N=300 because it seemed to find independent 
sequences with scores higher than 25 the most often:

              Independent sequences
Cycle length  (score > 20) (score > 25)  Complete Cycles  Num Steps
============  =========================  ===============  =========
         200             6            1               10       2000
         300             5            3                6       2000
         400             4            1                5       2000
         500             4            1                4       2000
         600             3            1                3       2000

Second, I chose a high temperature.  I tried T=2,3,4,5,6 and I scaled the cycle 
length to keep the slop constant.  T=2,3 didn't look very good at escaping 
maxima, and this was born out by auto-correlation times that were longer than 
the cycles themselves.  Between the remaining temperatures, I chose to move 
forward with T=5 for the same reason as above: it found independent sequences 
with scores above 25 the most often.

             Independent sequences
Temperature  (score > 20) (score > 25)  Complete Cycles  Num Steps
===========  =========================  ===============  =========
          4             4            1                8       2000
          5             5            3                6       2000
          6             4            1                5       2000

I'm aware that these results might be a fluke.  I used the same random seed for 
all these simulations, so the N=300, T=5 simulation was really only run once.  
It did slightly better than its peers in that run, but with <10 independent 
samples it's not really possible to say much.  That said, I don't think the 
precise parameters are crucial.  I just need to be in a region of parameter 
space where things work.

