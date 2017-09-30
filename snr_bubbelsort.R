
# bubbelsort snr lists ----------------------------------------------------

larger = function(pair) {
  print(pair)
  if(pair$snr[1] > pair$snr[2]) return(TRUE) else return(FALSE)
}

swap_if_larger = function(pair) {
  if(larger(pair)) {
    return(rev(pair)) 
  } else {
    return(pair)
  }
}

swap_pass = function(vec) { 
  for(i in seq(1, length(vec$snr)-1)) {
    vec[i:(i+1)] = swap_if_larger(vec[i:(i+1)])
  }
  return(vec)
}

bubble_sort = function(vec) {
  new_vec = swap_pass(vec)
  if(isTRUE(all.equal(vec, new_vec))) { 
    return(new_vec) 
  } else {
    return(bubble_sort(new_vec))
  }
}
