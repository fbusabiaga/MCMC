module read_input
"Read input file."

function read_input_file(name)

  # Define comment symbols
  comment_symbol = "#"

  # Create empty dictionary
  # options = Dict{String,Any}
  options = Dict()

  # Open file
  f_handle = open(name, "r")

  while ! eof(f_handle)
    line = readline(f_handle)
    n = findfirst(comment_symbol, line)
    
    # Remove comments
    if !isnothing(n)
      line = line[1:n[1]-1]   
    end
    
    # Fill dictionary
    if !isempty(line)
      values = split(line, " ", keepempty=false)
      options[String(values[1])] = values[2:end]
    end
    
  end
  close(f_handle)

  return options
end

end
