This tool uses the points in a blendshape with at least one inbetween shape to extrapolate a smooth trajectory for each vertex. This is done using a quadratic fit in the case of one in between, since that provides 3 points.
For a higher number of inbetweens natural cubic fit is used. Pretty niche but one use case might be when you need to slide something over a curved surface, such as a blinking eyelid.

Supports blend weights outside the 0 to 1 range



https://github.com/maccollo/Maya-smooth-blendshape-inbetweens/assets/23036010/3fa01b81-468c-40c4-ab4a-c5a63b0f4583

