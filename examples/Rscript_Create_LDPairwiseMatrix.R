##########################################################################
# Script: Creating a LD pairwise SNP matrix (i.e. an input of varmodel3) #
# Author: Frederic Labbe                                                 #
# Date: 27/10/2021                                                       #
# Affiliation: University of Chicago                                     #
##########################################################################


############################################################
# Define the function that create a LD pairwise SNP matrix #
############################################################

ldmat <- function(nb_snps, nb_groups, nb_linked_snps_group, ld_type_groups) {
  unliked_snps = 1:nb_snps                                       # Create a list of unlinked SNPs.
  pairwise_ld = matrix(0, nb_snps, nb_snps)                      # Create a matrix full of zeros with a number of rows and columns equal to the number of SNPs.
  diag(pairwise_ld) = 1.0                                        # Replace the diagonal values of the matrix by 1 (i.e. max LD).
  nb_linked_pair_snps = 0.0                                      # Define the number of paired linked SNPs as 0.
  if (min(nb_linked_snps_group) > 1) {                           # Test if the number of linked SNPs in each group is higher than 1.
    if (sum(nb_linked_snps_group) <= nb_snps) {                  # Test if the number of linked SNPs in the groups is not too high.
      if (length(ld_type_groups) == nb_groups) {                 # Test if the number of group LD intensity match the number of groups.
        for (i in 1:nb_groups) {                                 # For each linked SNP group:
          type = ld_type_groups[i]                               # Define the intensity of LD in this group.
          nb_linked_snps = nb_linked_snps_group[i]               # Define the number of linked SNPs.
          linked_snps = sample(unliked_snps, nb_linked_snps)     # Sample X SNPs from the unlinked SNPs, which will become the linked SNPs.
          linked_snps = sort(linked_snps)                        # Sort the linked SNPs (from lowest to highest).
          unliked_snps = setdiff(unliked_snps, linked_snps)      # Update the list of unliked SNPs.
          j = linked_snps[1]                                     # Define j as the first (i.e. smallest) linked SNP in the group (e.g. SNP 2 in group [2, 4, 6]).
          for (k in linked_snps) {                               # For each linked SNP:
            if (j == k) {                                        # If the SNP k and SNP j are identical:
              pairwise_ld[k, j] = as.double(1.0)                 # Define the pairwise LD between SNP k and SNP j as equals to 1 (e.g. LD between SNP 2 and SNP 2 is 1.0).
            } else {                                             # Otherwise (i.e.):
              nb_linked_pair_snps = nb_linked_pair_snps + 1      # Increment the number of paired linked SNPs.
              if (type == "high") {                              # If high intensity of LD in this SNP group:
                pairwise_ld[k, j] = round(runif(1, 0.7, 1.0), 2) # Define the pairwise high LD between SNP k and SNP j.
              } else if (type == "moderate") {                   # If moderate intensity of LD in this SNP group:
                pairwise_ld[k, j] = round(runif(1, 0.4, 0.7), 2) # Define the pairwise moderate LD between SNP k and SNP j.
              } else if (type == "low") {                        # If low intensity of LD in this SNP group:
                pairwise_ld[k, j] = round(runif(1, 0.1, 0.4), 2) # Define the pairwise low LD between SNP k and SNP j.
              } else {                                           # Otherwise (i.e. undefined LD intensity):
                print("Error: unknown intensity of LD!")         # Return an error.
              }
            }
          }
        }
        pairwise_ld[upper.tri(pairwise_ld)] <- t(pairwise_ld)[upper.tri(pairwise_ld)]
      } else {
        print("Error: the number of group LD intensity does not match the number of groups!")
      }
    } else {
      print("Error: the number of linked SNPs in the groups is too high (i.e. higher than the total number of SNPs!")
    }
  } else {
    print("Error: the number of linked SNPs in a group should be higher than 1!")
  }
  output = list(pairwise_ld, nb_linked_pair_snps)
  return(output)
}


#################################
# Test and define the arguments #
#################################

# 1st example of command line: Rscript --vanilla Rscript_Create_LDPairwiseMatrix.R 24 2 'c(2, 4)' 'c("high", "moderate")'
# 2nd example of command line: Rscript --vanilla Rscript_Create_LDPairwiseMatrix.R 24 2 'c(2)' 'c("low")' 1
args = commandArgs(trailingOnly = TRUE) # This function scans the arguments which have been supplied when the current R session was invoked. 
if (length(args) < 4) {                 # Test if there is four argument: if not, return an error.
  stop("the first four arguments must at least be supplied:\n
       1. Number of SNPs (required).\n
       2. Number of linked SNP group(s) (required).\n
       3. Number of linked SNPs in each group (required).\n
       4. Intensity of LD in each linked SNP group, i.e. low, moderate, or high (required).\n
       5. Suffix of the output file (optional).", call. = FALSE)
} 
nb_snps = as.numeric(args[1])                      # Number of SNPs.
nb_groups = as.numeric(args[2])                    # Number of linked SNP group(s).
nb_linked_snps_group = eval(parse(text = args[3])) # Number of linked SNPs in each group.
ld_type_groups = eval(parse(text = args[4]))       # Intensity of LD in each linked SNP group.
                                                   # Options: high (i.e. [0.7-1.0]), moderate (i.e. [0.4-0.7]), or low (i.e. [0.1-0.4]).
if (length(args) == 5) {                           # If the number of arguments is equals to 5:
  suffix = paste("_v", args[5], sep = "")          # Define the suffix of the output file as indicated in the arguments.
} else {                                           # Otherwise:
  suffix = ""                                      # Define the suffix of the output file as the default one.
}


#####################################
# Creating a LD pairwise SNP matrix #
#####################################

pairwise_ld = ldmat(nb_snps, nb_groups, nb_linked_snps_group, ld_type_groups) # Create the LD pairwise SNP matrix.
write.table(format(x = pairwise_ld[[1]], digits = 2),
            file = paste("pairwise_ld_", nb_snps, "snps_", nb_groups, "groups_", pairwise_ld[[2]], "pairs", suffix,".txt", sep = ""),
            row.names = F, col.names = F, quote = F, sep = "\t") # Export the LD pairwise SNP matrix.

