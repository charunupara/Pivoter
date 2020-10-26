using DelimitedFiles
using StatsBase

"""
pivoter.jl
Charun Upara

This file contains code that implements the algorithm from The Power of Pivoting for Exact Clique Counting by Shweta Jain and C. Seshadhri

The code is strongly modeled after the authors' implementation in C, found at https://bitbucket.org/sjain12/pivoter/

To run the code, call Julia from the command line with 3 arguments
    1. path to input file, stored as a .edges file
    2. the number of nodes in the graph
    3. the highest number of k to count k-cliques (use 0 to count cliques of all sizes)

    For example:
        julia pivoter.jl data/email-Enron.edges 36692 0

The code also requires the file 'nCr.txt' in data/ to function properly. 
The file simply stores pre-computed results from the binomial() function to speed up the overall runtime.

Lastly, because most data files are too large to attach over the internet, they can be found online by searching 
their name. To generate a .edges file, follow the instruction found at https://bitbucket.org/sjain12/pivoter/

Table of contents 
    1. Definitions of LinkedList data structures and their procedures (used for degeneracy ordering)
    2. Definitions of NeighborList data structures and code for computing degeneracy ordering of a graph
    3. Misc functions
    4. Helper functions for SCT building and clique listing
    5. Main functions
    6. Program entry point
"""

#### Linked List definition start####

# a node in a linked list, named link to follow the authors' naming conventions
mutable struct Link
    data::Union{Int64, Nothing} # allow for null value
    next::Union{Link, Nothing} # allow for null nodes
    prev::Union{Link, Nothing} # allow for null nodes
end

Link() = Link(nothing, nothing, nothing)

# linked list data structure
mutable struct LinkedList
    head::Union{Link, Nothing} 
    tail::Union{Link, Nothing} 
    length::Int64
end

LinkedList() = LinkedList(Link(), Link(), 0) # empty linked list


"""
Procedure: 
    * is_tail
Parameters: 
    * list, a Link
Purpose:
    * checks whether the link is the tail of a linked list
Produces:
    * result, a boolean
Preconditions:
    * list is not null
Postconditions:
    * result == true if list is tail. Otherwise, result == false
"""
function is_tail(list::Link)
    @assert list !== nothing
    return list.next === nothing
end


"""
Procedure: 
    * add_after
Parameters: 
    * list, a Link
    * data, an integer
Purpose:
    * add a new link after the current link
Produces:
    * new_link, a Link
Preconditions:
    * list is not null
    * list.next is not null
Postconditions:
    * [as seen in code]
"""
function add_after(list::Link, data::Int64)
    @assert list !== Link()
    @assert list.next !== Link()
    new_link = Link(data, list.next, list)
    list.next.prev = new_link
    list.next = new_link

    return new_link
    
end


"""
Procedure: 
    * add_before
Parameters: 
    * list, a Link
    * data, an integer
Purpose:
    * add a new link before the current link
Produces:
    * new_link, a Link
Preconditions:
    * list is not null
    * list.prev is not null
Postconditions:
    * [as seen in code]
"""
function add_before(list::Link, data::Int64)
    @assert list !== Link()
    @assert list.prev !== Link()
    new_link = Link(data, list, list.prev)
    list.prev.next = new_link
    list.prev = new_link

    return new_link
end

"""
Procedure: 
    * remove_link
Parameters: 
    * list, a Link
Purpose:
    * removes references to list from other Links
Produces:
    * list, a Link
Preconditions:
    * list is not null
    * list.next is not null
    * list.prev is not null
Postconditions:
    * [as seen in code]
"""
function remove_link(list::Link)
    @assert list !== Link()
    @assert list.next !== Link()
    @assert list.prev !== Link()

    list.next.prev = list.prev
    list.prev.next = list.next
    
    list.next = nothing
    list.prev = nothing

    return list
end

"""
Procedure: 
    * delete_link
Parameters: 
    * list, a Link
Purpose:
    * checks preconditions and gets data before calling remove_link()
Produces:
    * data, a integer
Preconditions:
    * list is not null
    * list.next is not null
    * list.prev is not null
Postconditions:
    * data = list.data
    * no Link is connected to list
"""
function delete_link(list::Link)
    @assert list !== Link()
    @assert list.next !== Link()
    @assert list.prev !== Link()

    data = list.data

    list = remove_link(list)

    return data

end



"""
Procedure: 
    * create_linked_list
Parameters: 
Purpose:
    * creates an empty linked list
Produces:
    * linked_list, a LinkedList
Preconditions:
Postconditions:
    * all fields is null or 0 except
        * linked_list.head.next = linked_list.tail
        * linked_list.tail.prev = linked_list.head
        to make it circular
"""
function create_linked_list()
    linked_list = LinkedList(Link(), Link(), 0)
    linked_list.head.prev = nothing
    linked_list.head.next = linked_list.tail
    linked_list.head.data = nothing

    linked_list.tail.prev = linked_list.head
    linked_list.tail.next = nothing
    linked_list.tail.data = nothing

    return linked_list
end

"""
Procedure: 
    * is_empty
Parameters: 
    * linked_list, a LinkedList
Purpose:
    * checks if a LinkedList is empty
Produces:
    * result, a boolean
Preconditions:
Postconditions:
    * result = true is linked list is empty. Otherwise, result = false
"""
function is_empty(linkedlist::LinkedList)
    return is_tail(linkedlist.head.next)
end

"""
Procedure: 
    * add_first
Parameters: 
    * linked_list, a LinkedList
    * data, an integer
Purpose:
    * add a Link to the beginning of linked_list
Produces:
    * result, a Link
Preconditions:
Postconditions:
    * linked_list.head.next = result
    * linked_list.length += 1
"""
function add_first(linkedlist::LinkedList, data::Int64)
    linkedlist.length += 1
    return add_after(linkedlist.head, data)
end

"""
Procedure: 
    * add_last
Parameters: 
    * linked_list, a LinkedList
    * data, an integer
Purpose:
    * add a Link to the end of linked_list
Produces:
    * result, a Link
Preconditions:
Postconditions:
    * linked_list.tail.prev = result
    * linked_list.length += 1
"""
function add_last(linkedlist::LinkedList, data::Int64)
    linkedlist.length += 1
    return add_before(linkedlist.tail, data)
end

"""
Procedure: 
    * get_first
Parameters: 
    * linked_list, a LinkedList
Purpose:
    * gets the first element of linked_list
Produces:
    * result, a Link
Preconditions:
    * linked_list is not empty
Postconditions:
"""
function get_first(linkedlist::LinkedList)
    @assert(!is_empty(linkedlist))
    return linkedlist.head.next.data
end

#### Linked List definition end ####










#### Degeneracy order definition start ####

# Data structure for each node in finding degeneracy orientation
mutable struct NeighborList
    vertex::Union{Int64, Nothing}
    earlier::Union{LinkedList, Nothing}
    later::Union{LinkedList, Nothing}
    order_number::Union{Int64, Nothing}
end

NeighborList() = NeighborList(nothing, nothing, nothing, nothing) # empty neighbor list

# An extension of NeighborList for the entire vertex set
mutable struct NeighborListArray
    vertex::Union{Int64, Nothing}
    earlier::Union{Vector{Int64}, Nothing}
    earlier_degree::Union{Int64, Nothing}
    later::Union{Vector{Int64}, Nothing}
    later_degree::Union{Int64, Nothing}
    order_number::Union{Int64, Nothing}
end

NeighborListArray() = NeighborListArray(nothing, nothing, nothing, nothing, nothing, nothing) # empty neighbor list array


"""
Procedure: 
    * compute_degeneracy_order_array
Parameters: 
    * list, a linked list
    * size, an integer
Purpose:
    * finds degeneracy orientation of a graph, the first step in building the SCT
Produces:
    * ordering_array, a linked list
Preconditions:
    * list represents the adjacency list of the input graph
    * length(list) == size
Postconditions:
    * ordering_array represents the degeneracy ordering of the input graph
"""
function compute_degeneracy_order_array(list::Vector{LinkedList}, size::Int64)
    # This procedure is not described in the pseudocode; code adapted from authors' implementation in C
    ordering = Vector{NeighborList}()

    vertices_by_degree = Vector{LinkedList}()
    vertex_locator = Vector{Link}()
    degree = Vector{Int64}()

    for i in 1:size
        push!(vertices_by_degree, create_linked_list())
        push!(ordering, NeighborList())
        ordering[i].earlier = create_linked_list()
        ordering[i].later = create_linked_list()
    end

    for i in 1:size
        # + 1 because index of the degree is used to access elements of other arrays, and Julia uses 1-indexing
        push!(degree, (list[i].length) + 1)
        to_add = add_first(vertices_by_degree[degree[i]], i)
        push!(vertex_locator, to_add)  
    end

    
    cur_degree = 1 # same reason as above
    num_vertices_removed = 0



    # compute ordering for each vertex
    while num_vertices_removed < size
        if !is_empty(vertices_by_degree[cur_degree])

            vertex = get_first(vertices_by_degree[cur_degree])
            delete_link(vertex_locator[vertex])
            ordering[vertex].vertex = vertex
            ordering[vertex].order_number = num_vertices_removed + 1
            degree[vertex] = -1

            neighbor_list = list[vertex]
            neighbor_link = neighbor_list.head.next

            while !is_tail(neighbor_link)
                neighbor = neighbor_link.data
                
                if degree[neighbor] != -1 
                    delete_link(vertex_locator[neighbor])
                    add_last(ordering[vertex].later, neighbor)

                    degree[neighbor] -= 1
                    if degree[neighbor] != -1
                        vertex_locator[neighbor] = add_first(vertices_by_degree[degree[neighbor]], neighbor)
                    end
                else
                    add_last(ordering[vertex].earlier, neighbor)
                end
                neighbor_link = neighbor_link.next
            end
            num_vertices_removed += 1
            cur_degree = 1
        else
            cur_degree += 1
        end
    end

    
    ordering_array = Vector{NeighborListArray}()


    for i in 1:size
        push!(ordering_array, NeighborListArray())
    end

    # moves results for each vertex to create the final ordering array
    for i in 1:size
        ordering_array[i].vertex = ordering[i].vertex
        ordering_array[i].order_number = ordering[i].order_number
        ordering_array[i].later_degree = ordering[i].later.length


        ordering_array[i].later = zeros(ordering_array[i].later_degree)

        j=1
        cur = ordering[i].later.head.next
        while !is_tail(cur)
            ordering_array[i].later[j] = cur.data
            cur = cur.next
            j+=1
        end

        ordering_array[i].earlier_degree = ordering[i].earlier.length
        ordering_array[i].earlier = zeros(ordering_array[i].earlier_degree)

        j=1
        cur = ordering[i].earlier.head.next

        while(!is_tail(cur))
            ordering_array[i].earlier[j] = cur.data
            cur = cur.next
            j += 1
        end
    end



    return ordering_array
end

#### Degeneracy order definition end ####







#### Misc functions start ####


"""
Procedure: 
    * read_graph
Parameters: 
    * filename, a string
    * n, an integer
Purpose:
    * parses input graph from a .edges file to an adjacency list
Produces:
    * adj_list, a linked list
Preconditions:
    * filename is a valid .edges file
    * n is the correct number of vertices for the corresponding input file
    * max(filename) + 1 == n
Postconditions:
    * adj_list represents the correct adjacency list for the input graph
"""
function read_graph(filename::AbstractString, n::Int64)
    adj_list = Vector{LinkedList}()
    for i in 1:n
        push!(adj_list, create_linked_list())
    end
    E = readdlm(filename,header=true)[1]
    E = Int.(E)'
    for i = 1:2:length(E)
        u = E[i] + 1
        v = E[i+1] + 1

        add_last(adj_list[u], v)
        add_last(adj_list[v], u)
    end

    return adj_list
end

"""
Procedure: 
    * populate_ncr
Parameters: 
Purpose:
    * reads in the nCr file generated by the authors that contains pre-generated values for binomial(drop, i)
        that will be used in incrementing clique counts
Produces:
    * ncr, a 1001x401 matrix 
Preconditions:
    * nCr.txt is a valid and correct file
Postconditions:
    * ncr stores the contents of nCr.txt
"""
function read_ncr()
    ncr = readdlm("data/nCr.txt", ',', Float64)
    return ncr
end

# used as a global variable to reduce the work of passing it with every function call
global ncr = nothing
ncr = read_ncr()



#### Misc functions end ####


#### Helper functions for SCT building and clique listing start ####


"""
Procedure: 
    * fill_p_and_x
Parameters: 
    * vertex, an integer
    * order_number, an integer
    * vertex_sets, an array of integers
    * vertex_lookup, an array of integers
    * ordering_array, an array of NeighborListArray
    * begin_x, an integer
    * begin_p, an integer
    * begin_r, an integer
    * new_begin_x, nothing
    * new_begin_p, nothing
    * new_begin_r, nothing
Purpose:
    * moves a vertex into R and define sets P and X as its earlier and later neighbors, respectively
        (sets P, R, and X are defined by their respective indices in vertex_sets;
         for example, set X = vertex_sets[begin_x] to vertex_sets[begin_p].)
Produces:
Preconditions:
    * [no other]
Postconditions:
    * vertex is at location begin_r
    * vertex set and associated arrays represent an ordering of vertices described above
    * new_begin_x, new_begin_p, new_begin_r are valid indices for the sets
"""
function fill_p_and_x(vertex, order_number, vertex_sets, vertex_lookup, ordering_array, neighbors_in_p, num_neighbors, begin_x, begin_p, begin_r, new_begin_x, new_begin_p, new_begin_r)
   vertex_location = vertex_lookup[vertex]
   begin_r -= 1 

   # move vertex to R
   vertex_sets[vertex_location] = vertex_sets[begin_r]
   vertex_lookup[vertex_sets[begin_r]] = vertex_location
   vertex_sets[begin_r] = vertex
   vertex_lookup[vertex] = begin_r
   new_begin_r = begin_r
   new_begin_p = begin_r 


   # find index where earlier neighbors start
   j = 1
   while j <= ordering_array[order_number].later_degree
    neighbor = ordering_array[order_number].later[j]
    neighbor_location = vertex_lookup[neighbor]
    new_begin_p -= 1 

    vertex_sets[neighbor_location] = vertex_sets[new_begin_p]
    vertex_lookup[vertex_sets[new_begin_p]] = neighbor_location
    vertex_sets[new_begin_p] = neighbor
    vertex_lookup[neighbor] = new_begin_p
    j+=1
   end



   new_begin_x = new_begin_p 

   # set the amount of possible neighbors of each vertex that are in P
   j = new_begin_p
   while j <= new_begin_r
    vertex_in_p = vertex_sets[j]
    num_neighbors[vertex_in_p] = 1
    arr_size = min((new_begin_r-new_begin_p)+1, (ordering_array[vertex_in_p].later_degree+ordering_array[vertex_in_p].earlier_degree))
    neighbors_in_p[vertex_in_p] = ones(arr_size)
    j+=1
   end


   # set neighbors in P and update its count for each vertex
   j = new_begin_p
   while j <= new_begin_r
    vertex_in_p = vertex_sets[j]
    k = 1
    while k <= ordering_array[vertex_in_p].later_degree
       later_neighbor = ordering_array[vertex_in_p].later[k] 
       later_neighbor_location = vertex_lookup[later_neighbor]

       
       if later_neighbor_location >= new_begin_p && later_neighbor_location < new_begin_r

            if num_neighbors[vertex_in_p] > length(neighbors_in_p[vertex_in_p])
                push!(neighbors_in_p[vertex_in_p], later_neighbor)
                num_neighbors[later_neighbor] +=1
            else
                neighbors_in_p[vertex_in_p][num_neighbors[vertex_in_p]] = later_neighbor
                num_neighbors[vertex_in_p] +=1
            end
           
            if num_neighbors[later_neighbor] > length(neighbors_in_p[later_neighbor])
                push!(neighbors_in_p[later_neighbor], vertex_in_p)
                num_neighbors[later_neighbor] +=1
            else
                neighbors_in_p[later_neighbor][num_neighbors[later_neighbor]] = vertex_in_p
        
                num_neighbors[later_neighbor] +=1     
            end
     
       end

       k+=1
    end

    j+=1
   end

   return begin_x, begin_p, begin_r,new_begin_x, new_begin_p, new_begin_r
end


"""
Procedure: 
    * find_best_pivot
Parameters: 
    * vertex_sets, an array of integers
    * vertex_lookup, an array of integers
    * neighbors_in_p, an array of arrays
    * num_neighbors, an array of integers
    * begin_x, an integer
    * begin_p, an integer
    * begin_r, an integer
Purpose:
    * finds pivot vertex for the process of building the SCT
Produces:
    * pivot, an integer
    * pivot_non_neighbors, an array of integers
    * num_non_neighbors, an integers
Preconditions:
    * [no other]
Postconditions:
    * pivot is the vertex with the largest neighborhood in P
    * pivot_non_neighbors contains non-neighbors of pivot
    * num_non_neighbors = length(pivot_non_neighbors)
"""
function find_best_pivot(vertex_sets, vertex_lookup, neighbors_in_p, num_neighbors, begin_x, begin_p, begin_r)
    pivot = 0
    max_intersect_size = 0


    j = begin_p

    # count number of neighbors in P for each vertex
    while j <= begin_r
        vertex = vertex_sets[j]
        num_potential_neighbors = min((begin_r-begin_p)+1, num_neighbors[vertex])
        num_neighbors_in_p = 0
        k = 1 
        while k < num_potential_neighbors
            neighbor = neighbors_in_p[vertex][k]
            neighbor_location = vertex_lookup[neighbor]

            if neighbor_location >= begin_p && neighbor_location < begin_r
                num_neighbors_in_p += 1
            else
                break
            end
            k += 1
        end

        # update pivot based on neighborhood size
        if num_neighbors_in_p >= max_intersect_size
            pivot = vertex
            max_intersect_size = num_neighbors_in_p
        end
        j += 1
    end


    pivot_non_neighbors = zeros((begin_r-begin_p) + 1)

    for i in 1:(begin_r-begin_p)+1
        pivot_non_neighbors[i] = vertex_sets[(begin_p-1)+i]
    end

   
    
    num_non_neighbors = (begin_r-begin_p)+1

    num_pivot_neighbors = min(num_non_neighbors, num_neighbors[pivot])


    # find and mark non-neighbors
    j = 1
    while j < num_pivot_neighbors
        neighbor = neighbors_in_p[pivot][j]
        neighbor_location = vertex_lookup[neighbor]
        
        if neighbor_location >= begin_p && neighbor_location < begin_r
            loc = neighbor_location-begin_p
            pivot_non_neighbors[loc+1] = -1
        else
            break
        end

        j += 1
        
    end
    j = 1
    while j < num_non_neighbors
        vertex = pivot_non_neighbors[j]
        if vertex == -1
            pivot_non_neighbors[j] = pivot_non_neighbors[num_non_neighbors]
            num_non_neighbors -= 1
            continue
        end
        j+=1
    end

    

    filter!(x->x≠-1, pivot_non_neighbors) # use to resolve a bug
    unique!(pivot_non_neighbors) # use to resolve a bug
    return pivot, pivot_non_neighbors,length(pivot_non_neighbors)
end



"""
Procedure: 
    * move_to_r
Parameters: 
    * vertex, an integer
    * vertex_sets, an array of integers
    * vertex_lookup, an array of integers
    * neighbors_in_p, an array of arrays
    * num_neighbors, an array of integers
    * begin_x, an integer
    * begin_p, an integer
    * begin_r, an integer
    * new_begin_x, an integer
    * new_begin_p, an integer
    * new_begin_r, an integer
Purpose:
    * moves a vertex to R and updates P and X
Produces:
    * new_begin_x, an integer
    * new_begin_p, an integer
    * new_begin_r, an integer
Preconditions:
    * [no other]
Postconditions:
    * new_begin_x, new_begin_p, new_begin_r represent new beginning indices of each set
    * vertex is moved to set R and sets P and X are adjusted accordingly
"""
function move_to_r(vertex, vertex_sets, vertex_lookup, neighbors_in_p, num_neighbors, begin_x, begin_p, begin_r)
    vertex_location = vertex_lookup[vertex]
    
    # move to R
    begin_r -= 1
    vertex_sets[vertex_location] = vertex_sets[begin_r]
    vertex_lookup[vertex_sets[begin_r]] = vertex_location
    vertex_sets[begin_r] = vertex
    vertex_lookup[vertex] = begin_r


    new_begin_x = begin_p
    new_begin_p = begin_p
    new_begin_r = begin_p
    

    size_of_p = (begin_r - begin_p)+1

    # define new locations for each set
    j = begin_p
    while j <= begin_r
        neighbor = vertex_sets[j]
        neighbor_location = j
        num_potential_neighbors = min(size_of_p, num_neighbors[neighbor])

        k = 1
        while k <= num_potential_neighbors
            if neighbors_in_p[neighbor][k] == vertex && new_begin_r <= begin_r
                vertex_sets[neighbor_location] = vertex_sets[new_begin_r]
                vertex_lookup[vertex_sets[new_begin_r]] = neighbor_location
                vertex_sets[new_begin_r] = neighbor
                vertex_lookup[neighbor] = new_begin_r
                new_begin_r += 1
            end
            k += 1
        end
        j += 1
    end

    # update each vertex's neighbors in P
    j = new_begin_p
    while j <= new_begin_r
        this_vertex = vertex_sets[j]
        num_potential_neighbors = min(size_of_p, num_neighbors[this_vertex])
        num_neighbors_in_p = 1

        k = 1

        while k <= num_potential_neighbors
            neighbor = neighbors_in_p[this_vertex][k]
            neighbor_location = vertex_lookup[neighbor]
            if neighbor_location >= new_begin_p && neighbor_location < new_begin_r
                neighbors_in_p[this_vertex][k] = neighbors_in_p[this_vertex][num_neighbors_in_p]
                neighbors_in_p[this_vertex][num_neighbors_in_p] = neighbor
                num_neighbors_in_p += 1
            end
            k += 1
        end
        j += 1
    end

    
    return begin_x, begin_p, begin_r, new_begin_x, new_begin_p, new_begin_r
end


"""
Procedure: 
    * move_r_to_x
Parameters: 
    * vertex, an integer
    * vertex_sets, an array of integers
    * vertex_lookup, an array of integers
    * begin_x, an integer
    * begin_p, an integer
    * begin_r, an integer
Purpose:
    * move a vertex from R to X
Produces:
Preconditions:
    * [no other]
Postconditions:
    * vertex is moved to set X and the indices are updated correctly
"""
function move_r_to_x(vertex, vertex_sets, vertex_lookup, begin_x, begin_p, begin_r)
    vertex_location = vertex_lookup[vertex]

    vertex_sets[vertex_location] = vertex_sets[begin_p]
    vertex_lookup[vertex_sets[begin_p]] = vertex_location
    vertex_sets[begin_p] = vertex
    vertex_lookup[vertex] = begin_p

    begin_p += 1
    begin_r += 1

    return begin_x, begin_p, begin_r
    
end


#### Helper functions for SCT building and clique listing end ####



#### Main functions start ###


"""
Procedure: 
    * list_cliques_recursive
Parameters: 
    * clique_counts, an array of Float64
    * vertex_sets, an array of integers
    * vertex_lookup, an array of integers
    * neighbors_in_p, an array of arrays
    * num_neighbors, an array of integers
    * begin_x, an integer
    * begin_p, an integer
    * begin_r, an integer
    * max_k, an integer
    * rsize, an integer
    * drop, an integer
Purpose:
    * finishes building the SCT and updates clique counts
Produces:
Preconditions:
    * vertex_sets is an array of length n
    * vertex_lookup contains the location of vertex_lookup[i]
    * neighbors_in_p may be an array of empty arrays or the neighbors of a node that are in the set P
    * num_neighbors represents the size of each node's neighborhood
    * begin_x, begin_p, and begin_r correspond to the initial index of each set in the vertex set 
Postconditions:
    * [no other]
"""
function list_cliques_recursive(clique_counts, vertex_sets, vertex_lookup, neighbors_in_p, num_neighbors, begin_x, begin_p, begin_r, max_k, rsize, drop)

    # condition is checking whether the SCT has been completely built
    # if built, update clique counts accordingly then exit the recursive call
    if ((begin_p >= begin_r) || (rsize-drop>max_k))
        i = drop
        while i >= 0 && (rsize-i <= max_k)
            k = rsize - i
            # nCr[i,j] is simply binomial(i,j) that is predefined in a seperate file to reduce the run time
            clique_counts[k+1] += ncr[drop+1,i+1]
            i -= 1 
        end
        return 
    end


    # if SCT is not built, we continue looking for the next pivot
    pivot, my_candidates_to_iterate,num_candidates_to_iterate = find_best_pivot(vertex_sets, vertex_lookup, neighbors_in_p, num_neighbors, begin_x, begin_p, begin_r)
    
    
    # iterate through pivot non-neighbors and recursively call the function on each
    if num_candidates_to_iterate > 0
        iterator = 1
        while iterator <= num_candidates_to_iterate
    
            vertex = my_candidates_to_iterate[iterator]
            vertex = convert(Int64, vertex)



            begin_x, begin_p, begin_r, new_begin_x, new_begin_p, new_begin_r = move_to_r(vertex, vertex_sets, vertex_lookup, neighbors_in_p, num_neighbors, begin_x, begin_p, begin_r)
            
            # recursive case
            if vertex == pivot
                list_cliques_recursive(clique_counts, vertex_sets, vertex_lookup, neighbors_in_p, num_neighbors, new_begin_x, new_begin_p, new_begin_r, max_k, rsize+1, drop+1)
            else
                list_cliques_recursive(clique_counts, vertex_sets, vertex_lookup, neighbors_in_p, num_neighbors, new_begin_x, new_begin_p, new_begin_r, max_k, rsize+1, drop)
            end
            

            begin_x, begin_p, begin_p = move_r_to_x(vertex, vertex_sets, vertex_lookup, begin_x, begin_p, begin_r)
        
            iterator +=1
        end

        # update P after recursive call
        iterator = 1
        while iterator <= num_candidates_to_iterate
            vertex = my_candidates_to_iterate[iterator]
            vertex = convert(Int64, vertex)
            vertex_location = vertex_lookup[vertex]

            begin_p -= 1
            vertex_sets[vertex_location] = vertex_sets[begin_p]
            vertex_sets[begin_p] = vertex
            vertex_lookup[vertex] = begin_p
            vertex_lookup[vertex_sets[vertex_location]] = vertex_location
            iterator += 1
        end
    end
    return 
end


"""
Procedure: 
    * list_cliques
Parameters: 
    * clique_counts, an array of Float64
    * ordering_array, an array of linked list
    * n, an integer
    * max_k, an integer
Purpose:
    * starts building SCT and calls recursive code to finish and update clique counts accordingly
Produces:
    * clique_counts, an array of Float64
Preconditions:
    * [no other]
Postconditions:
    * clique_counts[i] represents the number of i-1 cliques in the graph
"""
function list_cliques(clique_counts, ordering_array, n, max_k)


    vertex_sets = Vector{Int64}()
    vertex_lookup = Vector{Int64}()

    
    neighbors_in_p = Vector{Vector{Int64}}()
    num_neighbors = Vector{Int64}()

    for i in 1:n
        push!(vertex_sets, i)
        push!(vertex_lookup, i)
        push!(neighbors_in_p, [])
        push!(num_neighbors, 1)
    end



    begin_x = 0
    begin_p = 0
    begin_r = n

    # for each vertex
    for i in 1:n
        vertex = ordering_array[i].vertex
        new_begin_x = nothing
        new_begin_p = nothing
        new_begin_r = nothing

        # define sets P and X for current vertex
        begin_x, begin_p, begin_r, new_begin_x, new_begin_p, new_begin_r = fill_p_and_x(i, vertex, vertex_sets, vertex_lookup, ordering_array, neighbors_in_p, num_neighbors, begin_x, begin_p, begin_r, new_begin_x, new_begin_p, new_begin_r)
        
        # variables used to increment clique count by binomial(drop, drop-i)
        drop = 0
        rsize = 1

        list_cliques_recursive(clique_counts, vertex_sets, vertex_lookup, neighbors_in_p, num_neighbors, new_begin_x, new_begin_p, new_begin_r, max_k, rsize, drop)
        begin_r += 1
    end

    clique_counts[1] = 1
    return clique_counts
end





"""
Procedure: 
    * run_clique
Parameters: 
    * adj_list, an array of linked lists
    * n, an integer
    * max_k, an integer
Purpose:
    * declares appropriate variable before calling the main clique counting code and storing the results
Produces:
    * clique_counts, an array of Float64
    * max_k, an integer
Preconditions:
    * [no other]
Postconditions:
    * clique_counts[i] represents the number of i-1 cliques in the graph
    * max_k is updated if its initial value is 0
"""
function run_clique(adj_list::Array{LinkedList, 1}, n::Int64, max_k::Int64)
    
    
    deg = 0

    
    
    ordering_array = compute_degeneracy_order_array(adj_list, n)

    # finds the largest possible clique
    for i in 1:n
        if deg < ordering_array[i].later_degree
            deg = ordering_array[i].later_degree
        end

    end

    # if user wants to find all cliques, update max_k
    if max_k == 0
        max_k = deg + 1
    end

    
    clique_counts = zeros(max_k+1)

    clique_counts = list_cliques(clique_counts, ordering_array, n, max_k) # main work
   
    return clique_counts, max_k
    
end


"""
Procedure: 
    * count_clique
Parameters: 
    * adj_list, an array of linked lists
    * n, an integer
    * max_k, an integer
Purpose:
    * times and prints the results of clique counting code
Produces:
    * result, the amount of time in seconds the program takes to execute
Preconditions:
    * adj_list is the correct represtation of the input graph
    * length(adj_list) = n
    * max_k >= 0 
Postconditions:
    * clique_counts[i] represents the number of i-1 cliques in the graph
    * total_cliques = sum(clique_counts)
"""
function count_cliques(adj_list::Array{LinkedList, 1}, n::Int64, max_k::Int64)
    total_cliques = 0
    result = @timed run_clique(adj_list, n, max_k) #all main work inside start here, which is where timer starts
    clique_counts, max_k = result.value

    # print results
    println("Execution completed in ", result.time, " seconds.")
    for i in 1:max_k
        if clique_counts[i] != 0
            println(i-1, "-cliques: ", clique_counts[i])
            total_cliques += clique_counts[i]
        end
    end

    println("Total number of cliques: ", total_cliques)

    return result.time
    
end

#### Main functions end ####


#### Program entry point ####

# basic error checking
if length(ARGS) != 3
    error("Please pass in a '.edges file', the number of nodes in that network, and the largest clique size to count (pass in 0 to count all cliques).")
end

# parse command line arguments
graph = convert(String, ARGS[1])
n = parse(Int64, ARGS[2])
max_k = parse(Int64, ARGS[3])

# build adjacency list
graph_list = read_graph(graph, n)

# run code
count_cliques(graph_list, n, max_k)

