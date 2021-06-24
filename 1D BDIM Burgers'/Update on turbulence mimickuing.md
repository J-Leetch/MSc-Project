Main points from this week:

-	The 1D Burgers’ problem is good from the point of view of:
o	its simplicity and cheap computation;
o	convective-diffusive effects similar to N-S;
o	has wall boundaries naturally.

-	However, the lessons learned from it must be transferrable to 3D turbulence. Therefore, there are the following problems with the current test case:
o	the discontinuous forcing in time means that there is no turbulent “energy balance” and the fake turbulence is decaying. We want sustained turbulence.
o	
