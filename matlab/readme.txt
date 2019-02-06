DESCRIPTION:
- MATLAB algorithm that generates helicoidal phase masks, displays them in
  a transmission/reflection SLM and acquires images of the obtained vortex
  intensity patterns of a digital coronagraphy system. It is also possible
  to simulate the free-space propagation of two input starts (gaussian
  beams) and see the output.

PARTS:
- Analysis: principal scripts
- Data: the inputs of the algorithm are the acquired vortex images
- Data -> Datalogdir: specific measurement folder 
- Output: processed images or plots
- Tools: functions used in the program

INPUTS:
- See "Paramters.m"

OUTPUTS:
- Plots on the PC or on the SLMs
- Processed vortex images (or plots)

NOTES:
- Units: a.u (arbitrary units) and cm for lengths, radians for angles and
  um for wavelengths
- Helicoidal phase masks are inside the Laguerre-Gauss beams category
- The variable mask is complex and is wrappped: mask = exp(i*mask)
- All wrappped phases are shown on [-pi,pi] 
- Always execute the program whenver you are exactly inside its folder

CONTRIBUTORS:
Grupo de Óptica y Fotónica - Universidad de Antioquia
Grupo de Óptica Aplicada - Universidad EAFIT
-Samuel Plazas Escudero
-Juan Jose Cadavid
-Rene Restrepo Gomez

Created on January 2018 and edited until June 2019.
