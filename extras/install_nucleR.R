#  Copyright 2019 Barcelona Supercomputing Center
#  
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  
#      http://www.apache.org/licenses/LICENSE-2.0
#  
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

  
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('dplyr')
BiocManager::install('IRanges')
BiocManager::install('GenomicRanges')
BiocManager::install('ShortRead',ask = FALSE)
BiocManager::install('doParallel')
BiocManager::install('ggplot2')
BiocManager::install('magrittr')
BiocManager::install("nucleR")
