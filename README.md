This tool uses the points in a blendshape with at least one inbetween shape to extrapolate a smooth trajectory for each vertex. This is done using a quadratic fit in the case of one in between, since that provides 3 points.
For a higher number of inbetweens natural cubic fit is used. One use case might be when you need to slide something over a curved surface, such as a blinking eyelid.

The code is not very optimised and it takes a bit of time to generate the extra inbetweens, but still faster and more convenient than doing it by hand.

https://github.com/maccollo/Maya-smooth-blendshape-inbetweens/assets/23036010/8edfc10c-0913-43ee-928f-9fecdd0d5979

