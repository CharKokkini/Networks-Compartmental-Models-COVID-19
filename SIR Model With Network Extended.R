#Copyright 2024 Charikleia Kokkini

#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at

#http://www.apache.org/licenses/LICENSE-2.0

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

#-------------------------------------------------------------------------------

#Status Key:
# 0 - susceptible
# 1 - infectious
# 2 - recovered

library(blockmodels)
library(igraph)

#Graph generation subroutines
erdos_renyi <-function(n) {
  #Function to create preferential attachment graph 
    #Takes number of nodes as arguement
  
  g1 <- erdos.renyi.game(n = n, p = 0.01, type = "gnp", directed = FALSE )
  
  return(g1)
}

stochastic_block <- function(n) {
  #Function to create graph based on the Stochastic Block Model 
    #Will represent the ideal lockdown scenario
  
  prob_edge_between_comms <- function(new_n, num_comm){
    #Function to generate the probability of edges between communities
    return((nodes_per_comm / new_n)/sample(1:3, 1, replace = TRUE))
    }
  
  num_comm <- sample(180:190, 1, replace = TRUE)
  #Number of communities - random between x:y
  
  nodes_per_comm <- as.integer(n/num_comm)
  #Number of nodes per community
  
  new_n <- nodes_per_comm * num_comm
  #Remove the number of nodes that is the remainder 
  #Could range from 1 to 9 which is small relative to the total size so is inconsequential
  
  mat <- matrix(integer(), nrow = num_comm, ncol = num_comm)
  
  for (i in 1:num_comm){
    #i - row number
    for (j in i:num_comm){
      #j- column number
      if (i != j){
        value <- prob_edge_between_comms(new_n, num_comm)
        mat[i, j] <- value
        mat[j, i] <- value
        #This is to ensure the matrix is symmetric
        #Which is necessary as the graph is not directed
      } else {
        mat[i, j] <- prob_edge_between_comms(new_n, num_comm) * sample(5:10, 1)
        #Probability of edges inside a community will be 5 to 10 times greater 
        #(number selected at random)
      }
    }
  }
  
  vector_nodes_per_block <- rep(nodes_per_comm, num_comm)
  
  g1 <- sample_sbm(new_n, pref.matrix = mat, block.sizes = vector_nodes_per_block , directed = FALSE, loops = FALSE)
  
  return(g1)
}

pref_attachment <- function(n) {
  #Function to create a graph based on the Barabasi-Albert model with preferential attachment
  
  g1 <- barabasi.game(n, m = 2, directed = FALSE, zero.appeal = 0.1)
    #n - number of nodes
    #m - number of edges each node has
    #zero.appeal - the likelihood a node forms an edge with a node of high edge density
  
  return(g1)
}

#Subroutine to generate the new row in the table
generate_row <- function(nodes, current_time_step){
  susceptible_nodes <- as.list(which(nodes == '0'))
  #Stores the indices of all susceptible nodes
  
  infectious_nodes <- as.list(which(nodes == '1'))
  #Stores the indices of all infectious nodes
  
  recovered_nodes <- as.list(which(nodes == '2'))
  #Stores the indices of all recovered nodes
  
  num_susceptible_nodes <- length(susceptible_nodes)
  num_infectious_nodes <- length(infectious_nodes)
  num_recovered_nodes <- length(recovered_nodes)
  #Get the lengths of each of the lists and store them in variables
  #Lengths of each of the variables are the number of nodes in each compartment
  
  new_row <- list(current_time_step, 
                  num_susceptible_nodes, 
                  num_infectious_nodes, 
                  num_recovered_nodes)
  #Create a list of current time step, number of susceptible, infectious and recovered nodes
  #New row to be added to the table
  
  return(new_row)
}


#Main code
n <- 750
total_time_steps <- 200
Time_remain_infectious <- as.integer(14)
Time_remain_recovered <- as.integer(90)
prob_infection <- 0.2
prob_recovery <- 0.3
#Arbitrarily chosen values

table <- data.frame("Time Step" = c("0"),
                    "Susceptible" = c(as.character(n-1)), 
                    "Infectious" = c("1"),
                    "Recovered" = c("0"))
#Create a table with columns 'Time Step', 'Susceptible', 'Infectious' and 'Recovered'
#Where a row will be added at every iteration


graph_type <- "S"

if (graph_type == "E"){
  g1 <- erdos_renyi(n)

} else if (graph_type == "P") {
  g1 <- pref_attachment(n)
  
} else if (graph_type == "S"){
  g1 <- stochastic_block(n)
  
} else {
  stop("Invalid input. Please run again")
}
#Code to generate desired type of graph


LO <- layout_nicely(g1)
#layout.nicely() function adjusts the layout of the graph


nodes <- list()
#List 'nodes' with the status of all nodes in 0, 1 or 2

for (i in 1:n){
  nodes <- append(nodes, "0")
}
# Set all nodes to 0 initially as they are susceptible

nodes <- replace(nodes, 1, "1")
#Replace the first node's status to 1 (infectious)

t_infectious <- vector("list", n)
#Create another list 't_infectious' to store how long each node has been infectious for
  #The list has n contents and all are NULL initially
t_infectious[] <- replicate(n, as.integer(0))
#Set all values within the list to 0

t_infectious <- replace(t_infectious, which(nodes == "1"), as.integer(0))
#Replace the values in the list 't_infectious' from NULL to 0 for all infectious nodes

t_recovered <- vector("list", n)
#Create another list 't_recovered' to store how long each node has been recovered for
  #The list has n contents and all are NULL initially

node_colors <- c("0" = "green", "1" = "red", "2" = "blue")
#Nodes of status 0 meaning susceptible are coloured in green
#Nodes of status 1 meaning infectious are coloured in red
#Nodes of status 2 meaning recovered are coloured in blue

plot(g1, vertex.color = node_colors[as.character(nodes)],
     main = paste("Social Network Based on the SIR Model - Time Step: 0"),
     vertex.size = 5,
     vertex.label = "", 
     layout = LO)
#Sys.sleep(1)  
#Plot initial graph

recovered_nodes <- integer(0)

for (current_time_step in 1:total_time_steps) {
  
  for (edge in E(g1)) {
    source_node <- ends(g1, edge)[1]
    target_node <- ends(g1, edge)[2]
    #Selects an edge in the graph and names each node as the 'source_node' or the 'target_node'
    
    if (nodes[source_node] == "1" & nodes[target_node] == "0" & source_node != target_node) {
      if (runif(1, 0, 5) < prob_infection) {
        #runif(1) generates a random number between 0 and 1
        nodes[target_node] <- "1"
        t_infectious[target_node] <- as.integer(0)
        #If the 'source_node' is infectious and the 'target_node' is susceptible,
          #then the target node has a chance to become infectious
        #If the generated number is less than the 'prob_infection' then the source node is infectious 
        #and its t_infectious (time infectious) is set to 0
        
        t_recovered[target_node] <- as.integer(Time_remain_infectious * -1)
        #Set the time the node has been recovered to -Time_remain_infectious so that as it increments,
          #it will reach 0 when the node is no longer infectious and then can keep incrementing 
      }
    }
    else if (nodes[source_node] == "0" & nodes[target_node] == "1" & source_node != target_node) {
      #Check if the target_node is infectious and the target_node is susceptible
        #If so, do the same thing as above
      
      if (runif(1, 0, 5) < prob_infection) {
        nodes[source_node] <- "1"
        t_infectious[source_node] <- as.integer(0)
        
        t_recovered[target_node] <- as.integer(Time_remain_infectious * -1)
        #Set the time the node has been recovered to -Time_remain_infectious so that as it increments,
        #it will reach 0 when the node is no longer infectious and then can keep incrementing 
      }
    }
  }

  for (i in 1:n){
    c_node <- t_infectious[[i]]
    #c_node should store the value of how many time steps the i-th node has been infected
    
    if (c_node >= Time_remain_infectious & nodes[[i]] == "1" & runif(1) > prob_recovery) {
      #If the node has been infected for longer than the time specified in 'Time_remain_infectious'
      nodes[i] <- "2"
      #Set status to recovered
    }
      
    if (c_node >= (Time_remain_infectious + Time_remain_recovered) & nodes[[i]] == "2"){
      #If the node became infected as many time steps as specified above, then it has been infected 
        #for 14 time steps and recovered for 100
      nodes[i] <- "0"
      #Therefore, its status is set back to susceptible
    }
  }

  new_row <- generate_row(nodes, current_time_step)
  #Use the function generate_row to create a row to add to the table based on current data
    
  table[nrow(table) + 1,] <- new_row
  #Add the contents of 'new_row' as a new row in the table
  
  t_infectious <- lapply(t_infectious, function(x) if (!is.null(x)) as.integer(x + 1) else NULL)
  #Increment the contents of the 't_infectious' list by one UNLESS they are NULL in which case do nothing
  
  t_recovered <- lapply(t_recovered, function(x) if (!is.null(x)) as.integer(x + 1) else NULL)
  #Increment the contents of the 't_recovered' list by one UNLESS they are NULL in which case do nothing

  plot(g1, vertex.color = node_colors[as.character(nodes)],
       main = paste("Social Network Based on the SIR Model - Time Step:", current_time_step),
       vertex.size = 5,
       vertex.label = "",
       layout = LO)
  #Plot of updated graph
  
  legend("topright", legend = c("Susceptible", "Infected", "Recovered"),
         fill = node_colors, bty = "n", cex = 0.6)
  #Add key to the graphs
}

plot(table[["Susceptible"]], type = "l",
     col = "green",
     ylim= range(0: n),
     xlim = range(0: total_time_steps),
     xlab = "Time",
     ylab = "Population")

lines(table[["Infectious"]], col = "red")
lines(table[["Recovered"]], col = "blue")
#After all graphs have been plotted, plot a graph with lines representing
  #the number of susceptible, infectious and recovered nodes

legend("right", legend = c("Susceptible", "Infectious", "Recovered"),
       col = c("green", "red", "blue"), lty = 1, bty = "n", cex = 0.6)
#Add a key
  #bty = "n" - removes the box so it is not blocking the plot
  #cex changes the size of the letters





#THINGS TO DO
# # at the end means the action is done

#Generate pref attachment #
#Generate random graph with low and high edge density #
#Print every 10 time steps - print at end only #
#Make 500 time steps #
#Make 500 nodes #
#Keep infection rate the same for all graphs #
#Subroutine to generate each graphs #
#Add key of colours to graphs and plot #
#Figure out how to spread the nodes out
#Figure out how to stop the graph from changing shape #

#GRAPH (with axes) IS WEIRD - POSSIBLE SOLUTIONS:
#Reintroduce infection later on
  #Have a small probability that a node becomes infected without one of its neighbours 
    #being infected to account for the possibility that our graph has a finite number of
    #edges and cannot model all relationships
#Add probabilites to recovering- becoming susceptible again
  #Make the slopes gentler

#POTENTIALLY DO
#Add vaccinated compartment - make vaccinated quicker to recover
#Add death rate, birth rate and diseased compartment


#Ideally the relationships between the nodes would remain fixed between graphs
  #but the population is large enough to make the effect negligible 


#Low edge density - edge density = 0.005, prob_infection = 0.05


print(table)

file_path <- "C:\\Users\\xarak\\OneDrive\\School Stuff\\A-Level\\Project\\R\\Result Graphs\\SBM\\SBM10.txt"
write.table(table, file_path, sep="\t", row.names=FALSE)
#Save the table as a txt file

no_edges <- ecount(g1)

edge_dens <- edge_density(g1)

cat("Number of edges: ", no_edges)
cat("Edge density: ", edge_dens)
