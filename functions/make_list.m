function [list_i_all, list_j_all] = ...
    make_list(list_i_all,list_j_all,list_i, list_j)

global pointer

list_i_all(pointer + 1 : pointer + length(list_i))  = list_i;
list_j_all(pointer + 1 : pointer + length(list_j))  = list_j;
pointer = pointer + length(list_i);