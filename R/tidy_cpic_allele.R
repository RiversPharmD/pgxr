#' Tidy CPIC Allele Definition
#'
#' @param df the dataframe created by reading in a CPIC allele file. These
#'           can be found at: https://cpicpgx.org/guidelines/ and selecting a
#'           specific guideline. The file is a .xslx with the name
#'           "GENE_allele_definition_table" and is manually converted to a
#'           csv file prior to loading into R
#'
#' @return df: the tidied dataframe with six columns:
#'             gene: the gene that this file pertains to
#'             allele: the *# call assigned to the rs_number
#'             functional_status: what this variant does to gene function
#'             variant: the variant number assigned arbitrarily, as some
#'                      variants don't hav rs_numbers attached to them
#'             base: the base change causing this variant
#'             rs_number: the rs_number assigned to this variant
#'
#' @export
#'
#' @examples
#' tidy_cpic_allele(CFTR_allele_definition_table)
#'
tidy_cpic_allele <- function(df) {

  # Extracting Gene Name -------------------------------------------------------
  # removes the need of the user to enter the gene name separately

  df_name_as_string <- deparse(substitute(df))
  ## converts the variable df to a string "df"

  Gene <- str_remove(df_name_as_string, "_allele_definition_table")
  ## extracts the gene name, does require common formatting

  df <-  Filter(function(x) !(all(is.na(x))), df) ## eliminates the NA columns
  ## that are added by CPIC

  colnames(df) <- c("allele", "functional_status",
                    paste0("var_", 1:(ncol(df)-2)))
  ## renames the df with the appropriate titles

  df <- df %>%
    rowid_to_column()

  df <- df %>%
    filter(rowid <(str_which(df$allele, "NOTES:")-1)) %>%
    select(-rowid) %>%
    filter(functional_status != "Allele Functional Status" | is.na(functional_status))
  ## The is.na is necessary because filter automatically drops na's

  rs_id <- df %>%
    select(-allele) %>%
    filter(functional_status == "rsID") %>%
    gather(key = variant,
           value = rs_number,
           -functional_status) %>%
    mutate(rs_number = if_else(is.na(rs_number), "not_given", rs_number)) %>%
    select(-functional_status)
  ## creates a two column df that store the var_number and rs_number


  df <- df %>%
    filter(!is.na(allele)) ## gets rid of metadata

  df <- df %>%
    gather(key = variant,  ## melts data into a long structure
           value = base,
           -allele,
           -functional_status) %>%
    filter(!is.na(base)) %>%  ## removes the NA that are generated with the melt
    mutate(functional_status = if_else(is.na(functional_status), "Unknown",
                                       functional_status))


  df <- df %>%
    left_join(rs_id) %>%
    arrange(allele) %>%
    mutate(gene = Gene) %>%
    select(gene, everything())

  return(df)

}
