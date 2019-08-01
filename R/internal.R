# .onAttatch
#
# Welcome message for Autotuner
#
# Inspired by https://github.com/Huang-lab/oppti/blob/master/R/internal.R
.onAttach <- function (libname, pkgname){
    k <- paste0(
        "     ____________________________________________________________\n",
        paste0("    / __  / / / /__  __/ _   /__",
               "  __/ / / /\\\     / / ____/ __   / \n"),
        paste0("   / / / / / / /  / / / / /",
               " /  / / / / / /  \\\   / / /___/ /_/  /  \n"),
        paste0("  / /_/ / / / /  / / / / /",
               " /  / / / / / / /\\\ \\\ / / ____/  _   \\\ \n"),
        " / __  / /_/ /  / / / /_/ /  / / / /_/ / /  \\\   /  /__/  / /  /\n",
        "/_/ /_/_____/  /_/ /_____/  /_/ /_____/_/    \\\ /_____/__/ /__/   v.",
        utils::packageVersion( "Autotuner"),
        "\n",
        "https://github.com/crmclean/autotuner\n")
    
    packageStartupMessage(k)
    
    
}

