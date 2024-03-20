import maya.cmds as cmds
from time import time
from maya.api import OpenMaya as om
### linear algebra code from https://integratedmlai.com/system-of-equations-solution/
class LinAlgMatrix:
    "Wrapper class for the functions from https://integratedmlai.com/system-of-equations-solution/"
    @staticmethod
    def print_matrix(Title, M):
        print(Title)
        for row in M:
            print([round(x,3)+0 for x in row])
    @staticmethod            
    def print_matrices(Action, Title1, M1, Title2, M2):
        print(Action)
        print(Title1, '\t'*int(len(M1)/2)+"\t"*len(M1), Title2)
        for i in range(len(M1)):
            row1 = ['{0:+7.3f}'.format(x) for x in M1[i]]
            row2 = ['{0:+7.3f}'.format(x) for x in M2[i]]
            print(row1,'\t', row2)
    @staticmethod        
    def zeros_matrix(rows, cols):
        A = []
        for i in range(rows):
            A.append([])
            for j in range(cols):
                A[-1].append(0.0)

        return A
    @staticmethod
    def identity_matrix(n):
        I = LinAlgMatrix.zeros_matrix(n, n)
        for i in range(n):
            I[i][i] = 1.0

        return I
    @staticmethod
    def matrix_transpose(M):
        if not M: return []
        return [[M[i][j] for i in range(len(M))] for j in range(len(M[0]))]
    @staticmethod
    def copy_matrix(M):
        rows = len(M)
        cols = len(M[0])

        MC = LinAlgMatrix.zeros_matrix(rows, cols)

        for i in range(rows):
            for j in range(cols):
                MC[i][j] = M[i][j]

        return MC
    @staticmethod
    def matrix_multiply(A,B):
        rowsA = len(A)
        colsA = len(A[0])

        rowsB = len(B)
        colsB = len(B[0])

        if colsA != rowsB:
            print('Number of A columns must equal number of B rows.')
            return

        C = LinAlgMatrix.zeros_matrix(rowsA, colsB)

        for i in range(rowsA):
            for j in range(colsB):
                total = 0
                for ii in range(colsA):
                    total += A[i][ii] * B[ii][j]
                C[i][j] = total

        return C
    @staticmethod
    def linear_solve(A, B):
        AM = LinAlgMatrix.copy_matrix(A)
        n = len(A)
        BM = LinAlgMatrix.copy_matrix(B)
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

class Spline:
    """Class that creates a spline from a set of points. For 3 data points a quadratic spline will be created. For 4 or more a natural cubic spline will be created."""
    def __init__(self, x,y):
        if len(x) != len(y):
            raise ValueError("x and y must have the same length")
        if len(x) < 3:
            raise ValueError("spline must have at least 3 elements")
        cubicSpline = len(x) > 3

        if cubicSpline:
            a,b,c,d = self.calculate_cubic_spline_coefficients_with_quadratic_derivative_end(x, y)
            self.a,self.b,self.c,self.d = a,b,c,d
            self.x = x
        else:
            x_in = [0, x[1]-x[0], x[2]-x[0]]
            a,b,c = self.calculate_quadratic_coefficients(x_in, y)
            self.a, self.b, self.c, self.d = [0], [a], [b], [c]
            self.x = [x[0], x[-1]]
    def calculate_quadratic_coefficients(self,x,y):
        """Calculate the coefficients of a quadratic function that passes through the points (w0,y0), (w1,y1) and (w2, y2)."""
        y0,y1, y2 = y[0], y[1], y[2]
        x0, x1, x2 = x[0], x[1], x[2]
        print(x0, x1, x2)
        a_val = (y2*(x0-x1) + y1*(x2-x0) + y0*(x1-x2))/ (x0-x2) / (x1-x2) / (x0-x1)
        b_val = (y2*(x1**2 - x0**2) + y1*( x0**2 - x2**2) + y0*(x2**2 -x1**2)) / (x0 - x2)/(x1 - x2)/(x0 - x1) 
        c_val = (y0*(-x2 + x1)*x1 + x0*(y1 - x1**2*y2 + x0*(-y1 + x1*y2)))/((-x2 + x1)*x1 + x0*(x2 + x0*(-x2 + x1) - x1**2))
        return a_val, b_val, c_val

    def evaluate_spline(self,x_points):
        """This function evaluates the spline at the given x points and returns a list of y values."""
        n = len(self.x) - 1
        result = []
        for point in x_points:
            for i in range(n):
                if self.x[i] <= point <= self.x[i+1]:
                    dx = point - self.x[i]
                    value = self.a[i]*dx**3 + self.b[i]*dx**2 + self.c[i]*dx + self.d[i] 
                    result.append(value)
                    break
        return result
        
    def calculate_cubic_spline_coefficients_with_quadratic_derivative_end(self,x, y):
        """Generates the coefficients for a natural cubic spline where the boundary condition is that the derivative at the end points is equal to the derivative of a quadratic function passing through the first 3 and last 3 points."""
        n = len(x) - 1
        h = [x[i+1] - x[i] for i in range(n)]
        print(x)
        print(y)
        # Calculate quadratic coefficients for boundary conditions
        a0, b0, _ = self.calculate_quadratic_coefficients([x[0], x[1], x[2]], [y[0], y[1], y[2]])
        an, bn, _ = self.calculate_quadratic_coefficients([x[n-2], x[n-1], x[n]], [y[n-2], y[n-1], y[n]])
        
        # Calculate the first derivatives at the endpoints based on quadratic coefficients
        f_prime_0 = 2 * a0 * x[0] + b0
        f_prime_n = 2 * an * x[n] + bn
        
        # Set up the system of equations
        A = LinAlgMatrix.zeros_matrix(n + 1, n + 1)
        d = LinAlgMatrix.zeros_matrix(n + 1, 1)
        
        # Boundary conditions for the first derivatives
        A[0][0], A[n][n] = 2, 2
        A[0][1], A[n][n-1] = 1, 1  # This sets up the relationship for the derivative
        d[0][0], d[n][0] = 3 * ((y[1] - y[0]) / h[0] - f_prime_0), 3 * (f_prime_n - (y[n] - y[n-1]) / h[n-1])
        
        # Fill A and d for the internal points
        for i in range(1, n):
            A[i][i-1] = h[i-1]
            A[i][i] = 2 * (h[i-1] + h[i])
            A[i][i+1] = h[i]
            d[i][0] = 3 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
        
        # Solve the system for b
        b = LinAlgMatrix.linear_solve(A, d)
        b = LinAlgMatrix.matrix_transpose(b)[0]
        
        # Calculate c and a coefficients based on b
        c = [0]*n
        a = [0]*n
        for i in range(n):
            c[i] = (y[i+1] - y[i])/h[i] - h[i]*(b[i+1] + 2*b[i])/3
            a[i] = (b[i+1] - b[i]) / (3*h[i])
            
        # The d coefficients are simply y[i]
        d_coeff = [y[i] for i in range(n)]
        
        return a, b, c, d_coeff

class Inbetweener:
    """Wrapper class for the inbetweener tool."""
    @staticmethod
    def lerp(v0, v1, t):
        """Linear interpolation between two values."""
        return (1 - t) * v0 + t * v1
    
    @staticmethod
    def create_inbetween_shapes(base_shape,vertexPositions, weights, num_inbetweens):
        if len(weights) < 3:
            print("2 points dont make no curve!")
            return None, None
        weights_to_create = [Inbetweener.lerp(min(weights), max(weights), i / float(num_inbetweens + 1)) for i in range(1, num_inbetweens + 1)]
        shapes = Inbetweener.create_inbetween_shapes_from_splines(base_shape, vertexPositions, weights, num_inbetweens,weights_to_create)
        return weights_to_create, shapes
    
    @staticmethod
    def create_inbetween_shapes_from_splines(base_shape, vertex_positions, weights, num_inbetweens,weights_to_create):
        """Uses the given vertex positions and weights to generate a spline for each vertex, and from that generate a new mesh shape at the given weights_to_create. Returns a list of new mesh shapes"""
        num_vertices = int(len(vertex_positions[0]) / 3)

        # Pre-calculate the cubic spline coefficients for all vertices and dimensions
        splines_x = []
        splines_y = []
        splines_z = []
        #make new splinme
        for vtx_id in range(num_vertices):
            # Extract positions for current vertex across all weights
            x_positions = [vertex_positions[wi][vtx_id * 3] for wi in range(len(weights))]
            y_positions = [vertex_positions[wi][vtx_id * 3 + 1] for wi in range(len(weights))]
            z_positions = [vertex_positions[wi][vtx_id * 3 + 2] for wi in range(len(weights))]

            # Calculate and store spline coefficients for x, y, z
            splines_x.append(Spline(weights, x_positions))
            splines_y.append(Spline(weights, y_positions))
            splines_z.append(Spline(weights, z_positions))

        shapes = []
        

        for weight in weights_to_create:
            new_shape = cmds.duplicate(base_shape, name=f'inbetween_{weight}')[0]
            new_positions = []

            for vtx_id in range(num_vertices):
                # Evaluate splines at the current in-between weight for x, y, z
                new_x_position = splines_x[vtx_id].evaluate_spline([weight])[0]
                new_y_position = splines_y[vtx_id].evaluate_spline([weight])[0]
                new_z_position = splines_z[vtx_id].evaluate_spline([weight])[0]

                # Store new position for later bulk update
                new_positions.append([new_x_position, new_y_position, new_z_position])

            # Bulk update vertex positions using the optimized approach
            Inbetweener.set_vertex_positions_with_openmaya(new_shape, new_positions)
                
            shapes.append(new_shape)

        return shapes

    @staticmethod
    def set_vertex_positions_with_openmaya(mesh_name, new_positions):
        """Uses open maya to set new vertex positions for the given mesh. This is faster than using the cmds module."""
        selection_list = om.MSelectionList()
        selection_list.add(mesh_name)
        dagPath = selection_list.getDagPath(0)

        mfnMesh = om.MFnMesh(dagPath)
        points = om.MPointArray()

        for pos in new_positions:
            points.append(om.MPoint(pos[0], pos[1], pos[2]))

        mfnMesh.setPoints(points)

    @staticmethod
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

    @staticmethod
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
    
    @staticmethod
    def find_inbetween_weights_from_target_name(blendshape_node, target_name):
        """Finds and returns the inbetween weights, inbetween items (inbetween index for the node) and target index of a blendshape target"""
        target_index = Inbetweener.get_blendshape_target_index(blendshape_node, target_name)
        if target_index is not None:
            inbetween_weights, inbetween_items,target_item = Inbetweener.find_inbetween_weights_from_target_index(blendshape_node, target_index)
            if inbetween_weights:
                return inbetween_weights, inbetween_items,target_item,target_index
        return None

    @staticmethod
    def apply_inbetweens_to_blendshape(blendshape_node, base_shape,base_shape_copy, target_name, num_inbetweens):
        """
        Automatically detects the inbetween weight, calculates new inbetweens,
        and applies them to the specified blendshape target.
        """
        inbetween_weights, inbetween_items,target_item,target_index = Inbetweener.find_inbetween_weights_from_target_name(blendshape_node, target_name)
        vertex_positions = []
        weights = []
        base_mesh_positions = cmds.xform(f'{base_shape_copy}.vtx[*]', q=True, t=True, ws=True)
        if inbetween_weights is not None:
            Added_base_shape = False
            for weight, item in zip(inbetween_weights, inbetween_items):
                if weight > 0 and not Added_base_shape:
                    vertex_positions.append(base_mesh_positions)
                    weights.append(0)
                    Added_base_shape = True
                positions = Inbetweener.extract_inbetween_positions(blendshape_node,base_shape, target_index, item) 
                positions = [coordinate for sublist in positions for coordinate in sublist[:3]]
                positions = [base_mesh_positions[i] + positions[i] for i in range(len(base_mesh_positions))]
                vertex_positions.append(positions)
                weights.append(weight)
            target_shape = f'{target_name}_target'
            # Calculate and create new inbetweens
            time_start = time()
            weights, shapes = Inbetweener.create_inbetween_shapes(base_shape_copy,vertex_positions, weights, num_inbetweens)
            time_end = time()
            print("Time to calculate inbetweens: ", time_end - time_start)
            Inbetweener.add_inbetween_to_blendshape(base_shape,blendshape_node, target_name, shapes, weights,target_index)
            cmds.delete(shapes)
            # Apply new inbetweens to the blendshape node
            # This is where you'd integrate the newly created inbetween shapes into the blendshape node
            # as shown in the previous steps, which involves duplicating shapes and setting them as inbetweens
            # at calculated weights.
        else:
            print(f"No inbetween found for {target_name} in {blendshape_node}.")

    @staticmethod
    def extract_inbetween_positions(blendshape_node,base_shape, target_index, inbetween_index):
        """Takes a blendshape node, base shape, target index and inbetween index and returns the position offset of each vertex that is stored in the blendshape node."""

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
    
    @staticmethod
    def find_inbetween_weights_from_target_index(blendshape_node, target_index):
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

    @staticmethod
    def add_inbetween_to_blendshape(base_shape,blendshape_node, target_name, new_shapes, weights,target_index):
        """
        Adds a new inbetween shape to the specified blendshape target at the given weight.
        """
        for new_shape, weight in zip(new_shapes, weights):
        # Add the new shape as an inbetween at the specified weight
            #round the weight to 3 decimal places
            rounded_weight = round(weight, 3)
            cmds.blendShape(blendshape_node, edit=True, target=(base_shape, target_index, new_shape, rounded_weight), inBetween=True)

    @staticmethod
    def get_blendshape_nodes():
        """
        Returns a list of blendshape nodes in the scene.
        """
        return cmds.ls(type='blendShape')

    @staticmethod
    def get_blendshape_targets(blendshape_node):
        """
        Returns a list of targets for the given blendshape node.
        """
        return cmds.aliasAttr(blendshape_node, query=True)[::2]  # Skip every other entry to get target names

    @staticmethod
    def apply_inbetweens(*args):
        """
        Wrapper function to apply inbetweens based on UI selections.
        """
        blendshape_node = cmds.optionMenu("blendshapeNodeMenu", query=True, value=True)
        target_name = cmds.optionMenu("targetMenu", query=True, value=True)
        num_inbetweens = int(cmds.textField("numInbetweensField", query=True, text=True))
        base_shape = Inbetweener.get_base_shape_from_blendshape(blendshape_node)
        base_shape_no_deformation = Inbetweener.duplicate_without_deformation(base_shape,blendshape_node)[0]
        #reset the transforms of the base shape
        cmds.setAttr(base_shape_no_deformation + '.translate', 0,0,0)
        cmds.setAttr(base_shape_no_deformation + '.rotate', 0,0,0)
        cmds.setAttr(base_shape_no_deformation + '.scale', 1,1,1)
        Inbetweener.apply_inbetweens_to_blendshape(blendshape_node,base_shape, base_shape_no_deformation, target_name, num_inbetweens)
        cmds.delete(base_shape_no_deformation)
        cmds.inViewMessage(amg=f'Inbetweens generated for {target_name}.', pos='midCenter', fade=True)

    @staticmethod
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
                    return Inbetweener.get_base_shape_from_blendshape(connection)

        else:
            cmds.warning("No base shape found for blendshape node: " + blendshape_node)
            return None
    @staticmethod
    def update_target_menu(*args):
        """
        Updates the target menu based on the selected blendshape node.
        """
        blendshape_node = cmds.optionMenu("blendshapeNodeMenu", query=True, value=True)
        targets = Inbetweener.get_blendshape_targets(blendshape_node)
        
        # Retrieve the existing menu items and delete them
        menu_items = cmds.optionMenu("targetMenu", query=True, itemListLong=True)
        if menu_items:
            for item in menu_items:
                cmds.deleteUI(item, menuItem=True)
        
        # Populate the menu with new targets
        for target in targets:
            cmds.menuItem(label=target, parent="targetMenu")
    @staticmethod
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
        cmds.optionMenu("blendshapeNodeMenu", changeCommand=Inbetweener.update_target_menu)
        for node in Inbetweener.get_blendshape_nodes():
            cmds.menuItem(label=node)
        
        # Target menu (populated based on selected blendshape node)
        cmds.text(label="Select Target:")
        cmds.optionMenu("targetMenu")
        
        # Number of inbetweens
        cmds.text(label="Number of Inbetweens:")
        cmds.textField("numInbetweensField")
        
        # Generate button
        cmds.button(label="Generate Inbetweens", command=Inbetweener.apply_inbetweens)
        
        cmds.showWindow()

if __name__ == '__main__':
    Inbetweener.create_ui()
