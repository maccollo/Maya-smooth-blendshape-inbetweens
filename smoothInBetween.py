import maya.cmds as cmds
from time import time
from maya.api import OpenMaya as om
### linear algebra code from https://integratedmlai.com/system-of-equations-solution/
def print_matrix(Title, M):
    print(Title)
    for row in M:
        print([round(x,3)+0 for x in row])
        
def print_matrices(Action, Title1, M1, Title2, M2):
    print(Action)
    print(Title1, '\t'*int(len(M1)/2)+"\t"*len(M1), Title2)
    for i in range(len(M1)):
        row1 = ['{0:+7.3f}'.format(x) for x in M1[i]]
        row2 = ['{0:+7.3f}'.format(x) for x in M2[i]]
        print(row1,'\t', row2)
        
def zeros_matrix(rows, cols):
    A = []
    for i in range(rows):
        A.append([])
        for j in range(cols):
            A[-1].append(0.0)

    return A
def identity_matrix(n):
    I = zeros_matrix(n, n)
    for i in range(n):
        I[i][i] = 1.0

    return I

def copy_matrix(M):
    rows = len(M)
    cols = len(M[0])

    MC = zeros_matrix(rows, cols)

    for i in range(rows):
        for j in range(cols):
            MC[i][j] = M[i][j]

    return MC

def matrix_multiply(A,B):
    rowsA = len(A)
    colsA = len(A[0])

    rowsB = len(B)
    colsB = len(B[0])

    if colsA != rowsB:
        print('Number of A columns must equal number of B rows.')
        return

    C = zeros_matrix(rowsA, colsB)

    for i in range(rowsA):
        for j in range(colsB):
            total = 0
            for ii in range(colsA):
                total += A[i][ii] * B[ii][j]
            C[i][j] = total

    return C
def linear_solver(A, B):
    AM = copy_matrix(A)
    n = len(A)
    BM = copy_matrix(B)
    m = len(B[0])
    indices = list(range(n)) # allow flexible row referencing ***
    for fd in range(n): # fd stands for focus diagonal
        fdScaler = 1.0 / AM[fd][fd]
        # FIRST: scale fd row with fd inverse. 
        for j in range(n): # Use j to indicate column looping.
            AM[fd][j] *= fdScaler
            if j < m:
                BM[fd][j] *= fdScaler
        
        # SECOND: operate on all rows except fd row.
        for i in indices[0:fd] + indices[fd+1:]: # skip fd row.
            crScaler = AM[i][fd] # cr stands for current row
            for j in range(n): # cr - crScaler * fdRow.
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
                if j < m:
                    BM[i][j] = BM[i][j] - crScaler * BM[fd][j]
    return BM



def calculate_natural_cubic_spline_coefficients(x, y):
    n = len(x) - 1
    h = [x[i+1] - x[i] for i in range(n)]
    
    # Construct the matrix A for the system Ac = d
    A = zeros_matrix(n+1, n+1)
    A[0][0], A[n][n] = 1, 1
    for i in range(1, n):
        A[i][i-1] = h[i-1]
        A[i][i] = 2 * (h[i-1] + h[i])
        A[i][i+1] = h[i]
    # Construct the vector d
    d = zeros_matrix(n+1, 1)
    for i in range(1, n):
        d[i][0] = 3*((y[i+1] - y[i])/h[i] - (y[i] - y[i-1])/h[i-1])
    
    # Solve the system Ac = d for c
    c = linear_solver(A, d)
    
    # Calculate b and d coefficients based on c
    b = [0]*n
    d_coeff = [0]*n
    for i in range(n):
        b[i] = (y[i+1] - y[i])/h[i] - h[i]*(c[i+1][0] + 2*c[i][0])/3
        d_coeff[i] = (c[i+1][0] - c[i][0]) / (3*h[i])
        
    # a coefficients are simply y[i]
    a = [y[i] for i in range(n)]
    
    return a, b, c, d_coeff

def evaluate_spline(x, a, b, c, d, x_points):
    # Function to evaluate the spline at given points
    n = len(x) - 1
    result = []
    for point in x_points:
        for i in range(n):
            if x[i] <= point <= x[i+1]:
                dx = point - x[i]
                value = a[i] + b[i]*dx + c[i][0]*dx**2 + d[i]*dx**3
                result.append(value)
                break
    return result
def calculate_quadratic_coefficients_oneDimention(y0, y1, y2, w0, w1,w2):
    """Calculate the coefficients of a quadratic function that passes through the points (w0,y0), (w1,y1) and (w2, y2)."""
    a_val = (y2*(w0-w1) + y1*(w2-w0) + y0*(w1-w2))/ (w0-w2) / (w1-w2) / (w0-w1)
    b_val = (y2*(w1**2 - w0**2) + y1*( w0**2 - w2**2) + y0*(w2**2 -w1**2)) / (w0 - w2)/(w1 - w2)/(w0 - w1) 
    c_val = (y0*(-w2 + w1)*w1 + w0*(y1 - w1**2*y2 + w0*(-y1 + w1*y2)))/((-w2 + w1)*w1 + w0*(w2 + w0*(-w2 + w1) - w1**2))
        
    return a_val, b_val, c_val
def calculate_natural_cubic_spline_coefficients_with_quadratic_end(x, y):
    n = len(x) - 1
    h = [x[i+1] - x[i] for i in range(n)]
    
    # Calculate quadratic coefficients for boundary conditions
    a0, _, _ = calculate_quadratic_coefficients_oneDimention(y[0], y[1], y[2], x[0], x[1], x[2])
    an, _, _ = calculate_quadratic_coefficients_oneDimention(y[n-2], y[n-1], y[n], x[n-2], x[n-1], x[n])
    
    # Set second derivatives at the boundaries based on quadratic coefficients
    #(for some reason 1.4*a gave a better result than 2*a)
    s_dd_0 = 1.4 * a0
    s_dd_n = 1.4 * an
    
    # Construct the matrix A for the system Ac = d
    A = zeros_matrix(n+1, n+1)
    A[0][0], A[n][n] = 1, 1
    for i in range(1, n):
        A[i][i-1] = h[i-1]
        A[i][i] = 2 * (h[i-1] + h[i])
        A[i][i+1] = h[i]
    
    # Construct the vector d with the modified boundary conditions
    d = zeros_matrix(n+1, 1)
    d[0][0] = s_dd_0
    d[n][0] = s_dd_n
    for i in range(1, n):
        d[i][0] = 3*((y[i+1] - y[i])/h[i] - (y[i] - y[i-1])/h[i-1])
    
    # Solve the system Ac = d for c
    c = linear_solver(A, d)
    
    # Calculate b and d coefficients based on c
    b = [0]*n
    d_coeff = [0]*n
    for i in range(n):
        b[i] = (y[i+1] - y[i])/h[i] - h[i]*(c[i+1][0] + 2*c[i][0])/3
        d_coeff[i] = (c[i+1][0] - c[i][0]) / (3*h[i])
        
    # a coefficients are simply y[i]
    a = [y[i] for i in range(n)]
    
    return a, b, c, d_coeff
def calculate_quadratic_coefficients(base_pos, inbetween_pos, target_pos, base_weight_value, inbetween_weight,end_weight):
    a = []
    b = []
    c = []  # c is directly the base position, as it's y0 in our equations
    w0 = base_weight_value
    w1 = inbetween_weight
    w2 = end_weight
    # Calculate a and b using the derived formulas
    for y0, y1, y2 in zip(base_pos, inbetween_pos, target_pos):
        a_val = (y2*(w0-w1) + y1*(w2-w0) + y0*(w1-w2))/ (w0-w2) / (w1-w2) / (w0-w1)
        b_val = (y2*(w1**2 - w0**2) + y1*( w0**2 - w2**2) + y0*(w2**2 -w1**2)) / (w0 - w2)/(w1 - w2)/(w0 - w1) 
        c_val = (y0*(-w2 + w1)*w1 + w0*(y1 - w1**2*y2 + w0*(-y1 + w1*y2)))/((-w2 + w1)*w1 + w0*(w2 + w0*(-w2 + w1) - w1**2))
        
        a.append(a_val)
        b.append(b_val)
        c.append(c_val)
    return a, b, c
def lerp(v0, v1, t):
    return (1 - t) * v0 + t * v1

def create_inbetween_shapes(base_shape,vertexPositions, weights, num_inbetweens):
    start_time = time()
    if len(weights) < 3:
        #STOP
        print("2 points dont make no curve!")
        return None, None
    print(weights)
    if len(weights) == 3:
        #Since we only have 3 pointss we can use a quadratic curve
        return create_inbetween_shapes_quadratic_optimized(base_shape, vertexPositions, weights, num_inbetweens)
    # ITs CUBIC TIME
    print("4 or more points make for a cubic curve!")
    return create_inbetween_shapes_natural_cubic_optimized(base_shape, vertexPositions, weights, num_inbetweens)
    
def create_inbetween_shapes_natural_cubic_optimized(base_shape, vertexPositions, weights, num_inbetweens):
    num_vertices = int(len(vertexPositions[0]) / 3)

    # Pre-calculate the cubic spline coefficients for all vertices and dimensions
    spline_coeffs_x = []
    spline_coeffs_y = []
    spline_coeffs_z = []

    for vtx_id in range(num_vertices):
        # Extract positions for current vertex across all weights
        x_positions = [vertexPositions[wi][vtx_id * 3] for wi in range(len(weights))]
        y_positions = [vertexPositions[wi][vtx_id * 3 + 1] for wi in range(len(weights))]
        z_positions = [vertexPositions[wi][vtx_id * 3 + 2] for wi in range(len(weights))]

        # Calculate and store spline coefficients for x, y, z
        spline_coeffs_x.append(calculate_natural_cubic_spline_coefficients_with_quadratic_end(weights, x_positions))
        spline_coeffs_y.append(calculate_natural_cubic_spline_coefficients_with_quadratic_end(weights, y_positions))
        spline_coeffs_z.append(calculate_natural_cubic_spline_coefficients_with_quadratic_end(weights, z_positions))

    shapes = []
    weights_for_inbetweens = [lerp(min(weights), max(weights), i / float(num_inbetweens + 1)) for i in range(1, num_inbetweens + 1)]

    for weight in weights_for_inbetweens:
        new_shape = cmds.duplicate(base_shape, name=f'inbetween_{weight}')[0]
        new_positions = []

        for vtx_id in range(num_vertices):
            # Evaluate splines at the current in-between weight for x, y, z
            new_x_position = evaluate_spline(weights, *spline_coeffs_x[vtx_id], [weight])[0]
            new_y_position = evaluate_spline(weights, *spline_coeffs_y[vtx_id], [weight])[0]
            new_z_position = evaluate_spline(weights, *spline_coeffs_z[vtx_id], [weight])[0]

            # Store new position for later bulk update
            new_positions.append([new_x_position, new_y_position, new_z_position])

        # Bulk update vertex positions using the optimized approach
        set_vertex_positions_with_openmaya(new_shape, new_positions)
            
        shapes.append(new_shape)

    return weights_for_inbetweens, shapes
def create_inbetween_shapes_natural_cubic(base_shape, vertexPositions, weights, num_inbetweens):
    # Assuming calculate_natural_cubic_spline_coefficients is already defined
    # The number of vertices
    num_vertices = int(len(vertexPositions[0]) / 3)

    # Initialize the list to hold the new shapes
    shapes = []
    weights_for_inbetweens = [lerp(min(weights), max(weights), i / float(num_inbetweens + 1)) for i in range(1, num_inbetweens + 1)]
    
    for weight in weights_for_inbetweens:
        new_shape = cmds.duplicate(base_shape, name=f'inbetween_{weight}')[0]
        new_positions = []
        for vtx_id in range(num_vertices):
            # Extract x, y, z coordinates for the current vertex across all weights
            x_positions = [vertexPositions[wi][vtx_id*3] for wi in range(len(weights))]
            
            y_positions = [vertexPositions[wi][vtx_id*3 + 1] for wi in range(len(weights))]
            z_positions = [vertexPositions[wi][vtx_id*3 + 2] for wi in range(len(weights))]
            # Calculate spline coefficients for x, y, z
            x_coeffs = calculate_natural_cubic_spline_coefficients_with_quadratic_end(weights, x_positions)
            y_coeffs = calculate_natural_cubic_spline_coefficients_with_quadratic_end(weights, y_positions)
            z_coeffs = calculate_natural_cubic_spline_coefficients_with_quadratic_end(weights, z_positions)
            
            # Evaluate spline at the current in-between weight
            new_x_position = evaluate_spline(weights, x_coeffs[0], x_coeffs[1], x_coeffs[2], x_coeffs[3], [weight])[0]
            new_y_position = evaluate_spline(weights, y_coeffs[0], y_coeffs[1], y_coeffs[2], y_coeffs[3], [weight])[0]
            new_z_position = evaluate_spline(weights, z_coeffs[0], z_coeffs[1], z_coeffs[2], z_coeffs[3], [weight])[0]
            new_positions.append([new_x_position, new_y_position, new_z_position])
            # Move vertex to the new position
        #cmds.xform(f'{new_shape}.vtx[{vtx_id}]', t=(new_x_position, new_y_position, new_z_position), ws=True)
        set_vertex_positions_with_openmaya(new_shape,new_positions)
            
        shapes.append(new_shape)

    return weights_for_inbetweens, shapes
def get_skip_list(base_positions, inbetween_positions, target_positions):
    #If all three positions have the same x y and z coordinates we can skip the calculation
    skip_list = []
    for base, inbetween, target in zip(base_positions, inbetween_positions, target_positions):
        if base[0] == inbetween[0] == target[0] and base[1] == inbetween[1] == target[1] and base[2] == inbetween[2] == target[2]:
            skip_list.append(True)
        else:
            skip_list.append(False)
    return skip_list
def set_vertex_positions_with_openmaya(mesh_name, new_positions):
    selectionList = om.MSelectionList()
    selectionList.add(mesh_name)
    dagPath = selectionList.getDagPath(0)

    mfnMesh = om.MFnMesh(dagPath)
    points = om.MPointArray()

    for pos in new_positions:
        points.append(om.MPoint(pos[0], pos[1], pos[2]))

    mfnMesh.setPoints(points)
def create_inbetween_shapes_quadratic_optimized(base_shape, vertexPositions, weights, num_inbetweens):
    #debug the time it takes to calculate the inbetweens
    base_positions = [vertexPositions[0][i:i+3] for i in range(0, len(vertexPositions[0]), 3)]
    inbetween_positions = [vertexPositions[1][i:i+3] for i in range(0, len(vertexPositions[1]), 3)]
    target_positions = [vertexPositions[2][i:i+3] for i in range(0, len(vertexPositions[2]), 3)]
    skip_list = get_skip_list(base_positions, inbetween_positions, target_positions)
    base_weight = weights[0]
    inbetween_weight = weights[1]
    end_weight = weights[2]
    #print("GOTHERE")
    #print(skip_list)
    # Pre-calculate coefficients for all vertices
    coeffs = [calculate_quadratic_coefficients(base_pos, inbetween_pos, target_pos, base_weight, inbetween_weight, end_weight) 
              for base_pos, inbetween_pos, target_pos in zip(base_positions, inbetween_positions, target_positions)]
    weightsInterpolated = [lerp(base_weight, end_weight, i / float(num_inbetweens + 1)) for i in range(1, num_inbetweens + 1)]
    shapes = []

    for weight in weightsInterpolated:
        new_shape = cmds.duplicate(base_shape, name=f'inbetween_{weight}')[0]
        new_positions = []
        for vtx_id, (a, b, c) in enumerate(coeffs):
            new_positions.append([a[j]*(weight**2) + b[j]*weight + c[j] for j in range(3)])
            
            # Move vertices to new positions
            #This is slow as fuck
        set_vertex_positions_with_openmaya(new_shape,new_positions)
        shapes.append(new_shape)
    return weightsInterpolated, shapes
def create_inbetween_shapes_quadratic(base_shape,vertexPositions, weights, num_inbetweens):
    # Retrieve vertex positions for base, inbetween, and target
    #the positions are the offsets from base so we need to sum them with the base coordinates to get the true positions
    # Split positions into sublists of 3 (x, y, z)
    base_positions = vertexPositions[0]
    base_positions = [base_positions[i:i+3] for i in range(0, len(base_positions), 3)]
    inbetween_positions = vertexPositions[1]
    inbetween_positions = [inbetween_positions[i:i+3] for i in range(0, len(inbetween_positions), 3)]
    target_positions = vertexPositions[2]
    target_positions = [target_positions[i:i+3] for i in range(0, len(target_positions), 3)]

    base_weight = weights[0]
    inbetween_weight = weights[1]
    end_weight = weights[2]

    weightsInterpolated = []
    shapes = []
    for i in range(1, num_inbetweens + 1):
        weight =lerp(base_weight, end_weight, i / float(num_inbetweens + 1))
        new_shape = cmds.duplicate(base_shape, name=f'inbetween_{i}')[0]
        
        for vtx_id in range(len(base_positions)):
            a, b, c = calculate_quadratic_coefficients(base_positions[vtx_id], inbetween_positions[vtx_id], target_positions[vtx_id],base_weight, inbetween_weight,end_weight)
            
            new_positions = [a[j]*(weight**2) + b[j]*weight + c[j] for j in range(3)]
            
            # Move vertices to new positions
            cmds.xform(f'{new_shape}.vtx[{vtx_id}]', t=new_positions, ws=True)
        weightsInterpolated.append(weight)
        shapes.append(new_shape)
    return weightsInterpolated, shapes
def duplicate_without_deformation(mesh_name,blendshape_node, name = None, unlock_channels = True):
    """Turns off all deformers on the mesh, duplicates it, then turns the deformers back on and returns the duplicated mesh. The unlock channels arguments will unlock the transforms"""
    if not name:
        name = mesh_name + '#'
    history = cmds.listHistory(mesh_name)

    nodes_to_turn_back_on = []
    nodes_envelope_value = []
    for node in history:
        if cmds.attributeQuery('envelope', node=node, exists=True):
            nodes_to_turn_back_on.append(node)
            nodes_envelope_value.append(cmds.getAttr(f'{node}.envelope'))
            cmds.setAttr(f'{node}.envelope', 0)
            if node == blendshape_node:
                break

    duplicated_mesh = cmds.duplicate(mesh_name, name = name)
    #loop over the nodes to turn back on and reset them to their original value
    for i, node in enumerate(nodes_to_turn_back_on):
        cmds.setAttr(f'{node}.envelope', nodes_envelope_value[i])
    if unlock_channels:
        cmds.setAttr(duplicated_mesh[0] + '.tx', lock=False)
        cmds.setAttr(duplicated_mesh[0] + '.ty', lock=False)
        cmds.setAttr(duplicated_mesh[0] + '.tz', lock=False)
        cmds.setAttr(duplicated_mesh[0] + '.rx', lock=False)
        cmds.setAttr(duplicated_mesh[0] + '.ry', lock=False)
        cmds.setAttr(duplicated_mesh[0] + '.rz', lock=False)
        cmds.setAttr(duplicated_mesh[0] + '.sx', lock=False)
        cmds.setAttr(duplicated_mesh[0] + '.sy', lock=False)
        cmds.setAttr(duplicated_mesh[0] + '.sz', lock=False)
    return duplicated_mesh
def get_blendshape_target_names(blendshape_node):
    targets = cmds.listAttr(blendshape_node + '.w', multi=True)
    targetNames = []

    for target in targets:
        targetNames.append(target)
    return targetNames

def get_blendShape_target_connections(blendshape_node):
    targets = cmds.listAttr(blendshape_node + '.w', multi=True)
    target_connection_in = []
    target_connection_out = []
    for target in targets:
        incoming = cmds.listConnections(blendshape_node + '.' + target, source=True, destination=False, plugs=True)
        outgoing = cmds.listConnections(blendshape_node + '.' + target, source=False, destination=True, plugs=True)
        target_connection_in.append(incoming)
        target_connection_out.append(outgoing)
    return target_connection_in, target_connection_out

def break_blendShape_target_connections(blendshape_node):   
    targets = cmds.listAttr(blendshape_node + '.w', multi=True)
    for target in targets:
        incoming = cmds.listConnections(blendshape_node + '.' + target, source=True, destination=False, plugs=True)
        outgoing = cmds.listConnections(blendshape_node + '.' + target, source=False, destination=True, plugs=True)
        if incoming:
            for conn in incoming:
                cmds.disconnectAttr(conn, blendshape_node + '.' + target)
        if outgoing:
            for conn in outgoing:
                cmds.disconnectAttr(blendshape_node + '.' + target, conn)
def zero_blendsShape_target_weights(blendshape_node):   
    targets = cmds.listAttr(blendshape_node + '.w', multi=True)
    for target in targets:
        cmds.setAttr(blendshape_node + '.' + target, 0)
def bake_blendshape_painted_weights(blendshape_node):
    target_names = get_blendshape_target_names(blendshape_node)
    target_conections_in, target_conections_out = get_blendShape_target_connections(blendshape_node)
    break_blendShape_target_connections(blendshape_node)
    zero_blendsShape_target_weights(blendshape_node)

def get_blendshape_target_index(blendshape_node, target_name):
    """
    Takes the name of a blendshape target and returns the index.
    """
    alias_list = cmds.aliasAttr(blendshape_node, q=True)
    target_index = None
    for i, alias in enumerate(alias_list):
        if alias == target_name:
            target_index = i // 2  # Each target has two entries in the alias list (weight and input)
            print("Target name found")
            return target_index
    print("Target name not found")
    return None

def find_inbetween_weights(blendshape_node, target_name):
    target_index = get_blendshape_target_index(blendshape_node, target_name)
    if target_index is not None:
        inbetween_weights, inbetween_items,target_item = find_inbetween_weights2(blendshape_node, target_index)
        if inbetween_weights:
            return inbetween_weights, inbetween_items,target_item,target_index
    return None

def apply_inbetweens_to_blendshape(blendshape_node, base_shape,base_shape_copy, target_name, num_inbetweens):
    """
    Automatically detects the inbetween weight, calculates new inbetweens,
    and applies them to the specified blendshape target.
    """
    inbetween_weights, inbetween_items,target_item,target_index = find_inbetween_weights(blendshape_node, target_name)
    vertexPositions = []
    weights = []
    BaseMeshPositions = cmds.xform(f'{base_shape_copy}.vtx[*]', q=True, t=True, ws=True)
    if inbetween_weights is not None:
        Added_base_shape = False
        for weight, item in zip(inbetween_weights, inbetween_items):
            if weight > 0 and not Added_base_shape:
                vertexPositions.append(BaseMeshPositions)
                weights.append(0)
                Added_base_shape = True
            Positions = extract_inbetween_positions(blendshape_node,base_shape, target_index, item) 
            Positions = [coordinate for sublist in Positions for coordinate in sublist[:3]]
            Positions = [BaseMeshPositions[i] + Positions[i] for i in range(len(BaseMeshPositions))]
            vertexPositions.append(Positions)
            weights.append(weight)
        target_shape = f'{target_name}_target'
        # Calculate and create new inbetweens
        time_start = time()
        weights, shapes = create_inbetween_shapes(base_shape_copy,vertexPositions, weights, num_inbetweens)
        time_end = time()
        print("Time to calculate inbetweens: ", time_end - time_start)
        add_inbetween_to_blendshape(base_shape,blendshape_node, target_name, shapes, weights,target_index)
        cmds.delete(shapes)
        # Apply new inbetweens to the blendshape node
        # This is where you'd integrate the newly created inbetween shapes into the blendshape node
        # as shown in the previous steps, which involves duplicating shapes and setting them as inbetweens
        # at calculated weights.
    else:
        print(f"No inbetween found for {target_name} in {blendshape_node}.")
def extract_inbetween_positions(blendshape_node,base_shape, target_index, inbetween_index):
    total_vertices = cmds.polyEvaluate(base_shape, vertex=True)
    # Retrieve the list of components that have been modified
    inbetween_components_attr = f'{blendshape_node}.inputTarget[0].inputTargetGroup[{target_index}].inputTargetItem[{inbetween_index}].inputComponentsTarget'
    modified_components = cmds.getAttr(inbetween_components_attr)
    
    # Retrieve the modified positions
    inbetween_positions_attr = f'{blendshape_node}.inputTarget[0].inputTargetGroup[{target_index}].inputTargetItem[{inbetween_index}].inputPointsTarget'
    modified_positions = cmds.getAttr(inbetween_positions_attr)

    # Initialize a full list of positions with zeros
    full_positions = [[0.0, 0.0, 0.0, 1.0] for _ in range(total_vertices)]
    
    # Function to expand the vertex indices
    def expand_indices(index_string):
        if ":" in index_string:
            start, end = map(int, index_string.split(":"))
            return list(range(start, end + 1))
        else:
            return [int(index_string)]
    
    # Update the full positions list with modified positions
    current_modified_index = 0
    for component in modified_components:
        indices = expand_indices(component.split("[")[-1].split("]")[0])
        for index in indices:
            full_positions[index] = list(modified_positions[current_modified_index])
            current_modified_index += 1
    return full_positions
def find_inbetween_weights2(blendshape_node, target_index):
    """
    Attempt to find the inbetween weights for a specific target of a blendshape node by checking each
    inputTargetItem for available inbetween weights.
    """
    inbetween_weights = []
    inbetween_items = []
    # Attempt to query all inputTargetItems for the targetIndex
    target_items = cmds.getAttr(f'{blendshape_node}.inputTarget[0].inputTargetGroup[{target_index}].inputTargetItem', multiIndices=True)
    
    if target_items:
        for item in target_items:
            weight = (item - 5000) / 1000.0  # Convert the internal weight to the actual weight
            if True or 0 < weight and weight < 1:
                inbetween_weights.append(weight)
                inbetween_items.append(item)
            if weight == 1:
                target_item = item
    
    return inbetween_weights, inbetween_items, target_item
def add_inbetween_to_blendshape(base_shape,blendshape_node, target_name, new_shapes, weights,target_index):
    """
    Adds a new inbetween shape to the specified blendshape target at the given weight.
    """
    for new_shape, weight in zip(new_shapes, weights):
    # Add the new shape as an inbetween at the specified weight
        #round the weight to 3 decimal places
        rounded_weight = round(weight, 3)
        cmds.blendShape(blendshape_node, edit=True, target=(base_shape, target_index, new_shape, rounded_weight), inBetween=True)

def get_blendshape_nodes():
    """
    Returns a list of blendshape nodes in the scene.
    """
    return cmds.ls(type='blendShape')

def get_blendshape_targets(blendshape_node):
    """
    Returns a list of targets for the given blendshape node.
    """
    return cmds.aliasAttr(blendshape_node, query=True)[::2]  # Skip every other entry to get target names

def apply_inbetweens(*args):
    """
    Wrapper function to apply inbetweens based on UI selections.
    """
    blendshape_node = cmds.optionMenu("blendshapeNodeMenu", query=True, value=True)
    target_name = cmds.optionMenu("targetMenu", query=True, value=True)
    num_inbetweens = int(cmds.textField("numInbetweensField", query=True, text=True))
    base_shape = get_base_shape_from_blendshape(blendshape_node)
    base_shape_no_deformation = duplicate_without_deformation(base_shape,blendshape_node)[0]
    #reset the transforms of the base shape
    cmds.setAttr(base_shape_no_deformation + '.translate', 0,0,0)
    cmds.setAttr(base_shape_no_deformation + '.rotate', 0,0,0)
    cmds.setAttr(base_shape_no_deformation + '.scale', 1,1,1)
    apply_inbetweens_to_blendshape(blendshape_node,base_shape, base_shape_no_deformation, target_name, num_inbetweens)
    cmds.delete(base_shape_no_deformation)
    cmds.inViewMessage(amg=f'Inbetweens generated for {target_name}.', pos='midCenter', fade=True)
def get_base_shape_from_blendshape(blendshape_node):
    """
    Retrieves the base shape connected to the specified blendshape node.
    """
    # Find the geometry connected as input to the blendshape node
    connections = cmds.listConnections(blendshape_node + '.outputGeometry', source=False, destination=True)
    if connections:
        #check if the connection is a shape node
        for connection in connections:
            if cmds.nodeType(connection) == 'transform':
                return connection
            else:
                #If not we continue to the next connection
                return get_base_shape_from_blendshape(connection)

    else:
        cmds.warning("No base shape found for blendshape node: " + blendshape_node)
        return None
def update_target_menu(*args):
    """
    Updates the target menu based on the selected blendshape node.
    """
    blendshape_node = cmds.optionMenu("blendshapeNodeMenu", query=True, value=True)
    targets = get_blendshape_targets(blendshape_node)
    
    # Retrieve the existing menu items and delete them
    menu_items = cmds.optionMenu("targetMenu", query=True, itemListLong=True)
    if menu_items:
        for item in menu_items:
            cmds.deleteUI(item, menuItem=True)
    
    # Populate the menu with new targets
    for target in targets:
        cmds.menuItem(label=target, parent="targetMenu")

def create_ui():
    """
    Creates the UI for generating inbetweens.
    """
    window_id = 'generateInbetweensUI'
    
    if cmds.window(window_id, exists=True):
        cmds.deleteUI(window_id)
    
    cmds.window(window_id, title="Generate Inbetweens", widthHeight=(300, 150))
    cmds.columnLayout(adjustableColumn=True)
    cmds.text(label="Select Blendshape Node:")
    
    # Blendshape node menu
    cmds.optionMenu("blendshapeNodeMenu", changeCommand=update_target_menu)
    for node in get_blendshape_nodes():
        cmds.menuItem(label=node)
    
    # Target menu (populated based on selected blendshape node)
    cmds.text(label="Select Target:")
    cmds.optionMenu("targetMenu")
    
    # Number of inbetweens
    cmds.text(label="Number of Inbetweens:")
    cmds.textField("numInbetweensField")
    
    # Generate button
    cmds.button(label="Generate Inbetweens", command=apply_inbetweens)
    
    cmds.showWindow()


# Example of use
if __name__ == '__main__':
    create_ui()
    #blendshape_node = 'blendShape2'
    #base_shape = 'polySurface63'
    #target_name = 'polySurface68'
    #num_inbetweens = 10  # Number of new inbetweens to create

    #apply_inbetweens_to_blendshape(blendshape_node, base_shape, target_name, num_inbetweens)