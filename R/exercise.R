
rm(list = ls())

library(magrittr)

# read the poem into R

poem <-
  "My love is like a red, red rose that's newly sprung in June:
My love is like the melody that's sweetly played in tune.
As fair art thou, my bonnie lass, so deep in love am I:
And I will love thee still, my dear, till a' the seas gang dry.
Till a' the seas gang dry, my dear, and the rocks melt wi' the sun :
And I will love thee still, my dear, while the sands o' life shall run.
And fare thee weel, my only love, and fare thee weel a while!
And I will come again, my love, thou' it were ten thousand mile."

# tokenkize text using the penn treebank tokenizer

library(tokenizers)

# create a charter vector of tokens
(poem_tokenized <- tokenize_ptb(poem, simplify = TRUE))


# an alternative tokenizer from the tokenizer package

(poem_tokens <- unlist(tokenize_words(poem)))

# how big is the corpus: N
(N <- length(poem_tokens))


# locate a node of interest
node <- "love"

# how often does the node appear
(R_1 <- sum(poem_tokens == node))

# get the positions of the node
(node_where <- which(poem_tokens == node))

# define span: xL/xR
w <- 2

# get the collocate positions
collocates_where <- c(as.vector(outer(node_where, 1:w, "+")), as.vector(outer(node_where, 1:w, "-")))
sort(collocates_where)

# ensure that no position is smaller than 1 or larger than the text size
collocates_where <- collocates_where[collocates_where >= 1 & collocates_where <= length(poem_tokens)]
sort(collocates_where)

# get the tokens to the position
all_collocates <- poem_tokens[collocates_where]
all_collocates

# table it
O_11_named_vec <- table(all_collocates)
O_11_named_vec
# the above is, technically speaking, not a table but a named vector

# THAT would be a table
as.data.frame(O_11_named_vec)

# get just the names
O_11_names <- names(O_11_named_vec)

# pre-allocate space for the results vector
C_1_vec <- numeric(length = length(O_11_named_vec))

# loop over the collocates to determine how often they appear in the corpis
for (i in 1:length(O_11_named_vec)) {
  C_1_vec[i] <- sum(poem_tokens == O_11_names[i])
}

C_1_vec


as.data.frame(O_11_named_vec) %>%
  cbind(C_1_vec)

# Now we know
# the size of the corpus N
# the window size w
# the collocations of all collocates: O_11
# the frequency of all collocates: C_1
# the frequency of the node: R_1


# So we can calculate the measures 

E_11_vec <- (R_1 * C_1_vec) / N
E_11_vec


(MU_vec <- O_11_named_vec / E_11_vec)

sort(MU_vec, decreasing = TRUE)
