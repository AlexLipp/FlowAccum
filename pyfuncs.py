"""
These pure python functions are not used in the main code, but are included 
for reference. They are not used because they are slow compared to the
equivalent functions in the cython module. 
"""

import numpy
from typing import List


def d8_to_receivers(arr: numpy.ndarray) -> numpy.ndarray:
    """
    Converts a D8 flow direction array to a receiver array. The receiver array
    contains the index of the node that receives the flow from the current cell.
    All boundary cells are set to have themselves as receivers (sink).
    Elsewise, the receiver is set to the index of the node that receives the flow.
    A D8 of 1 means the right cell, 2 means the lower right, 4 means the bottom,
    eight means the lower left, 16 means the left, 32 means the upper left,
    64 means the top, and 128 means the upper right according to ESRI convention.
    """
    nrows = arr.shape[0]
    ncols = arr.shape[1]
    receivers = numpy.empty(arr.shape, dtype=int)
    receivers[:] = -1
    for i in range(nrows):
        for j in range(ncols):
            # Check if boundary cell
            if i == 0 or j == 0 or i == nrows - 1 or j == ncols - 1 or arr[i, j] == 0:
                receivers[i, j] = i * ncols + j
            elif arr[i, j] == 1:  # Right
                receivers[i, j] = i * ncols + j + 1
            elif arr[i, j] == 2:  # Lower right
                receivers[i, j] = (i + 1) * ncols + j + 1
            elif arr[i, j] == 4:  # Bottom
                receivers[i, j] = (i + 1) * ncols + j
            elif arr[i, j] == 8:  # Lower left
                receivers[i, j] = (i + 1) * ncols + j - 1
            elif arr[i, j] == 16:  # Left
                receivers[i, j] = i * ncols + j - 1
            elif arr[i, j] == 32:  # Upper left
                receivers[i, j] = (i - 1) * ncols + j - 1
            elif arr[i, j] == 64:  # Top
                receivers[i, j] = (i - 1) * ncols + j
            elif arr[i, j] == 128:  # Upper right
                receivers[i, j] = (i - 1) * ncols + j + 1
            else:
                raise ValueError(f"Invalid flow direction value: {arr[i, j]}")
    return receivers.flatten()


def count_donors(r: numpy.ndarray) -> numpy.ndarray:
    """
    Counts the number of donors that each cell has.

    Args:
        r: The receiver indices.

    """
    np = len(r)  # np = number of pixels
    d = numpy.zeros(np, dtype=int)
    for j in range(np):
        d[r[j]] += 1
    return d


def ndonors_to_delta(nd: numpy.ndarray) -> numpy.ndarray:
    """
    Converts a number of donors array to an index array that contains the location of where the list of
    donors to node i is stored.

    Args:
        nd: The "number of donors" array.

    """
    np = len(nd)  # np = number of pixels
    # Initialize the index array to the number of pixels
    delta = numpy.zeros(np + 1, dtype=int)
    delta.fill(np)
    delta[-2::-1] -= numpy.cumsum(nd[::-1])
    return delta


def ndonors_to_delta_slow(nd: numpy.ndarray) -> numpy.ndarray:
    """
    Converts a number of donors array to an index array that contains the location of where the list of
    donors to node i is stored.

    Args:
        nd: The donor array.

    """
    np = len(nd)
    # Initialize the index array to the number of pixels
    delta = numpy.zeros(np + 1, dtype=int)
    delta[np] = np
    for i in range(np, -1, -1):
        if i == np:
            continue
        delta[i] = delta[i + 1] - nd[i]

    return delta


def make_donor_array(r: numpy.ndarray, delta: numpy.ndarray) -> None:
    """
    Makes the array of donors. This is indexed according to the delta
    array. i.e., the donors to node i are stored in the range delta[i] to delta[i+1].
    So, to extract the donors to node i, you would do:
    donors[delta[i]:delta[i+1]]

    Args:
        r: The receiver indices.
        delta: The delta index array.

    """
    np = len(r)  # np = number of pixels
    # Define an integer working array w intialised to 0.
    w = numpy.zeros(np, dtype=int)
    # Donor array D
    D = numpy.zeros(np, dtype=int)
    for i in range(np):
        D[delta[r[i]] + w[r[i]]] = i
        w[r[i]] += 1

    return D


def add_to_stack(
    l: int, j: int, s: numpy.ndarray, delta: numpy.ndarray, donors: numpy.ndarray
):
    """
    Adds node l, and its donors (recursively), to the stack.

    Args:
        l: The node index.
        j: The stack index.
        s: The stack.
        delta: The index array.
        donor_info: The donor information array.

    """
    s[j] = l
    j += 1
    delta_l = delta[l]
    delta_lplus1 = delta[l + 1]

    for n in range(delta_l, delta_lplus1):
        m = donors[n]
        if m != l:
            j = add_to_stack(m, j, s, delta, donors)

    return j


def add_to_stack(
    l: int, j: int, s: numpy.ndarray, delta: numpy.ndarray, donors: numpy.ndarray
):
    """
    Adds node l, and its donors (recursively), to the stack.

    Args:
        l: The node index.
        j: The stack index.
        s: The stack.
        delta: The index array.
        donor_info: The donor information array.

    """
    s[j] = l
    j += 1
    delta_l = delta[l]
    delta_lplus1 = delta[l + 1]

    for n in range(delta_l, delta_lplus1):
        m = donors[n]
        if m != l:
            j = add_to_stack(m, j, s, delta, donors)

    return j


def make_donors_list_of_lists(
    D: numpy.ndarray, delta: numpy.ndarray
) -> List[List[int]]:
    """
    Makes the donor array as a list of lists. This is indexed according to the delta
    array. i.e., the donors to node i are stored in the range delta[i] to delta[i+1].
    So, to extract the donors to node i, you would do:
    donors[delta[i]:delta[i+1]]

    Args:
        D: The donor array.
        delta: The delta index array.

    """
    np = len(delta) - 1
    donors = [[] for _ in range(np)]
    for i in range(np):
        for j in range(delta[i], delta[i + 1]):
            donors[i].append(D[j])

    return donors


def build_stack(receivers: numpy.ndarray) -> numpy.ndarray:
    """
    Builds the stack of nodes in topological order, given the receiver array.
    Starts at the baselevel nodes and works upstream.

    Args:
        receivers: The receiver array (i.e., receiver[i] is the ID
        of the node that receives the flow from the i'th node).
    """
    orig = numpy.arange(len(receivers))
    baselevel_nodes = numpy.where(orig == receivers)[0]
    n_donors = count_donors(receivers)
    delta = ndonors_to_delta(n_donors)
    donors = make_donor_array(receivers, delta)
    stack = numpy.zeros(len(orig), dtype=int) - 1
    j = 0
    for b in baselevel_nodes:
        j = add_to_stack(b, j, stack, delta, donors)

    return stack


def accumulate_flow(
    receivers: numpy.ndarray, stack: numpy.ndarray, weights: numpy.ndarray
) -> numpy.ndarray:
    """
    Accumulates flow along the stack of nodes in topological order, given the receiver array,
    the ordered stack, and a weights array which contains the contribution from each node.

    Args:
        receivers: The receiver array (i.e., receiver[i] is the ID
        of the node that receives the flow from the i'th node).
        stack: The ordered stack of nodes.
        weights: The weights array (i.e., the contribution from each node).
    """
    n = len(receivers)
    accum = weights
    # Accumulate flow along the stack from upstream to downstream
    for i in range(n - 1, -1, -1):
        donor = stack[i]
        recvr = receivers[donor]
        if donor != recvr:
            accum[recvr] += accum[donor]

    return accum


def get_profile_segments(
    starting_nodes: numpy.ndarray,
    delta: numpy.ndarray,
    donors: numpy.ndarray,
    field: numpy.ndarray,
    threshold: float = 0,
) -> List[List[int]]:
    """
    Returns the channel segments for a D8 network, where each segment is a list of node indices. Only adds nodes to segments if some
    specified `field` value (e.g., upstream area) is greater or equal than the threshold. Each segment in the list of segments starts
    at a node and ends at a bifurcation or dead end. Each segment contains first the start node and then all nodes upstream of it
    in the order they are visited. The segments are ordered topologically, so that the first segment in the list
    is base-level, and the last segment in the list is an upstream tributary. Base level nodes present an edge case, and as such
    are always present *twice* in the list of segments. This prevents returning of segments containing *only* a single baselevel nodes.
    i.e., ensuring that every segment is a valid line (not a point).

    This function uses a stack (first in, first out) to keep track of nodes to visit. It also uses a stack to keep track of segments
    that are being built. This avoids recursion, which is slow.

    Args:
        starting_nodes: array of nodes to start from. Expects these to exceed threshold!
        delta: array of delta values
        donors: array of donor nodes
        field: array of field values
        threshold: threshold value for field

    Returns:
        List of segments, where each segment is a list of node indices
    """
    segments = []  # List of segments to return
    s = []  # FIFO Stack of nodes to visit
    seg_stack = []  # FIFO Stack of segments
    for b in starting_nodes:
        s.append(b)
        seg_stack.append([b])
    if len(s) == 0:
        # If there are no valid baselevel nodes, return an empty list
        return []
    curr_seg = seg_stack.pop()
    while len(s) > 0:
        node = s.pop()
        curr_seg.append(node)
        n_donors = 0
        # Loop over donors
        for n in range(delta[node], delta[node + 1]):
            m = donors[n]
            if m != node:
                # We don't want to add the node to the queue if it's the same as the current node
                if field[m] >= threshold:
                    # Only add the node to the queue if field > threshold.
                    s.append(m)
                    n_donors += 1
        if n_donors == 1:
            # We're still on the same segment, so we just continue...
            pass
        elif n_donors > 1:
            # We've reached a bifurcation! Add the current segment to the list of segments.
            segments.append(curr_seg)
            # Now we start a new segment for each donor, and put them in the segments queue.
            for _ in range(n_donors):
                seg_stack.append([node])
            # Pop the last element of the segment stack and continue from where we left off.
            curr_seg = seg_stack.pop()
        elif n_donors == 0:
            # We've reached a dead end! Add the current segment to the list of segments.
            segments.append(curr_seg)
            if len(seg_stack) == 0:
                # If the segments queue is empty, we're done!
                break
            else:
                # Otherwise, pop the last element of the segment stack and continue from where we left off.
                curr_seg = seg_stack.pop()
    return segments
