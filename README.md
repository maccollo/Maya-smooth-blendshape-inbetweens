This tool uses the points in a blendshape with at least one inbetween shape to extrapolate a smooth trajectory for each vertex. This is done by computing a smooth polynomial fit to the path of each vertex and then generating inbetween shapes from that path. Pretty niche but might be useful in some situations.

Supports blend weights outside the 0 to 1 range.

How to use: Download the script and run it. The tool window will pop up. Select the blendshape node, the blendshape target, and the number of inbetweens to generate.


Example 1: Making an eyelid move smoothly over an eyeball



https://github.com/maccollo/Maya-smooth-blendshape-inbetweens/assets/23036010/b2ca71ba-2d0c-4ffd-a9fd-d7d52c488e96



Example 2: Smoothing out the linear "kink" from a blendshape with a negative version of the shape 




https://github.com/maccollo/Maya-smooth-blendshape-inbetweens/assets/23036010/ed5c064d-7fe7-4266-b673-49e97d27288b

