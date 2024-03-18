This tool uses the points in a blendshape with at least one inbetween shape to extrapolate a smooth trajectory for each vertex. This is done by computing a smooth polynomial fit to the path of each vertex and then generating inbetween shapes from that path. Pretty niche but might be useful in some situations.

Supports blend weights outside the 0 to 1 range


Example 1: Making an eyelid move smoothly over an eyeball
https://github.com/maccollo/Maya-smooth-blendshape-inbetweens/assets/23036010/3fa01b81-468c-40c4-ab4a-c5a63b0f4583


Example 2: Smoothing out the linear "kink" from a blendshape with a negative version of the shape 
https://github.com/maccollo/Maya-smooth-blendshape-inbetweens/assets/23036010/296dd7fe-c673-45da-9394-f3853652f39a

