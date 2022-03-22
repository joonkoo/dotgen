# dotgen
MATLAB program for generating dot array stimuli for magnitude perception studies

'dotGenJP' and 'dotGenJP_peri' are MATLAB functions that creates magnitude dimension parameters for constructing dot arrays systematically across a 3D space defined by number (N), size (Sz), and spacing (Sp), according to the framework developed by DeWind, Adams, Platt, & Brannon (2015). The two dotGenJP functions use a slightly different definition of the Sz dimension, that is elaborated in Park, DeWind, Woldorff, & Brannon (2016). 

'dotField2GKA' is a MATLAB function that generates a field of non-overlapping dots of various radii, written by Geoffrey K. Adams. 

'script_generate_dots.m' illustrates how these two functions are used for a systematic construction of dot arrays. 

The code may contain errors, in which case please report them to the author. Please use them at your own risk.

For a reference to the theoretical aspect of this code, please cite DeWind et al. (2015) and Park et al. (2016). For a reference to the implementation, please cite Park (in press).

DeWind, N. K., Adams, G. K., Platt, M. L., & Brannon, E. M. (2015). Modeling the approximate number system to quantify the contribution of visual stimulus features. Cognition, 142, 247–265. https://doi.org/10.1016/j.cognition.2015.05.016

Park, J., DeWind, N. K., Woldorff, M. G., & Brannon, E. M. (2016). Rapid and Direct Encoding of Numerosity in the Visual Stream. Cerebral Cortex, 26(2), 748–763. https://doi.org/10.1093/cercor/bhv017

Park, J. (in press). Flawed stimulus design in additive-area heuristic studies. Cognition. https://doi.org/10.1016/j.cognition.2021.104919
