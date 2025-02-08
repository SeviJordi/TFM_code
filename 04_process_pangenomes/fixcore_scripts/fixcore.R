#!/usr/bin/env Rscript

# Script para curar alineamientos de genes
# Usage: Rscript eval_msa.R [Árbol] [vcf] [passing samples file] [filtering out samples file] [alineamiento]
args = commandArgs(trailingOnly=TRUE)

# LIBRERIAS ####
library("tidyr")
library("dplyr")
library("ape")
library("phylobase")
library("stringr")
library("phangorn")

# Inicio el dataframe removedf
removedf= data.frame(isolate=character(0))


# Obtener longitud del alineamiento, está en el vcf
msalen = grep("length", readLines(args[2]), value = T) %>% # longitud del alineamiento
  str_remove(., ".*=") %>% str_remove(., ">") %>% as.integer()

# leo el vcf de aa
vcf <- 
  read.csv(args[2], sep="\t", skip = 3,
           colClasses=c("REF"="character", "ALT"="character")) %>%  # para que no lea aa F y T como boleanos  
  select(-c(X.CHROM,QUAL,ID,FILTER,INFO,FORMAT)) %>%
  # cambio los numeros por el aa que tiene cada muestra
  mutate( across( .cols = c(ALT),  ~str_remove_all( ., ","))) %>%
  mutate( across( .cols = -c(POS, REF, ALT),  ~str_sub(ALT, .,.) ))  %>% 
  mutate( across( .cols = -c(POS, REF, ALT),  ~str_replace_all( ., "^$", REF ))) %>%
  select(-c(REF)) %>% # cambio la referencia que crea anp-sites por la emboss
  rename(REF = "EMBOSS_001")


vcf_list <- # paso el vcf de matriz a tabla MOLTEN
  vcf %>%
  select(-ALT) %>%
  pivot_longer(cols = -c( POS, REF), names_to = "isolate", values_to = "ALT") 

# Primero elimino todas aquellas muestras que tengan un 15% de gaps

tempo <- 
  vcf_list %>%  
  #Eliminamos las posiciones comunes con la REF
  filter(!(REF == "X" & ALT == "X")) %>%
  filter(ALT == "X") %>% 
  group_by(isolate) %>% 
  count() %>%
  filter(n > (0.10 * msalen)) %>% 
  select(isolate)

removedf = rbind(removedf, tempo)

testnvcf <- # miro si hay alguna muestra con m?s de 3 aa variantes seguidas
  vcf_list %>% 
  filter(REF != ALT) %>% 
  group_by(isolate) %>%
  mutate(GROUP = paste(isolate,cumsum(c(1L, diff(POS) > 1)), sep="-")) %>%
  group_by(isolate, GROUP) %>%
  #aQUI NECESITO TODA UNA EXPLICACION
  # diff(POS) TE LA LA DIFERENCIA ENTRE LA POS X Y X-1 (PORQUE EL LAG PREDETERMINADO ES 1, SE PUEDE CAmBIAR)
  # diff(POS) > 1 DEVUELVE TRUE/FALSE SI ENTRE UNA POSICION Y OTRA HAY MÁS DE 1
  # el 1L DE c(1L, diff(POS) > 1), CONVIERTE EL T/F A BINARIO
  # Y CUMSUM SUMA CUMULATIVAMENTE Y COMO ES BONARIO SE QUEDAN GRUPOS CONSECUTIVOS.
  filter(any(REF != "X") & any(ALT != "X") ) %>% # elimino todos los grupos que son gaps #TODO: preguntar a neris por esto
  #group_by(GROUP, ALT) %>%
  count() %>%
  filter(n >= 3 )  %>%
  group_by(isolate) %>%
  distinct(isolate)

if (nrow(testnvcf) >= 1){ # si hay m?s de 3 Pos variantes seguidas en el alineamiento
  
  # Leer el alineamiento
  fastaf <- read.FASTA(args[5])
  
 
  # Hacer un nj y quitamos ya los aislados con más de 15% de gaps  
  tree <- dist.dna(fastaf, pairwise.deletion = T) %>%
	replace(is.na(.), 0) %>%
        replace(is.infinite(.), 1) %>%
	nj() %>%
	drop.tip(pull(removedf, isolate))
  # Guardar el árbol en el fichero indicado
  #  write.tree(tree, file=args[1])
 

  L <- list(seq(0, tree$Nnode)) # Create node names

  tr <- makeNodeLabel(tree, "n", nodeList = L)  # Give names to nodes to convert to phylo4
  print("arbol guardado")
  # Las distancias que sean < 0 se sustituyen por 0
  tr$edge.length[tr$edge.length < 0] <- 0
  
  tr4df <- as(tr, "phylo4") %>% # convert tree to phylo4 S4 class
    as("data.frame") %>% # Extract dataframe
    filter(!is.na(edge.length)) # Eliminar la raiz
  print("distancias leidas")

  # Desviación dtandard para ramas terminales 
  sdmax_tip <- tr4df %>%  # extract 3 SD from terminal branches
    filter(node.type == "tip") %>% 
    summarise(sd = mean(edge.length)+(3*sd(edge.length))) %>% 
    pull(sd)
  print("desv est ram ext")
  # Desviación dtandard para ramas internas
  sdmax_int <- tr4df %>% # extract 3 SD from internal branches
    filter(node.type == "internal") %>% 
    summarise(sd = mean(edge.length)+(3*sd(edge.length))) %>%
    pull(sd)
  print("desv est ram int")
  # Búsqueda de outliers
  outliers_tips <- tr4df  %>% # select teminal outliers --> those nodes and tips above the 3SD 
    filter((node.type == "tip" & edge.length > sdmax_tip)) %>%
    mutate(Filter = "Test", Outlier = node)
  
  outliers_int <- tr4df  %>% # select internal outliers --> those nodes and tips above the 3SD 
    filter((node.type == "internal" & edge.length > sdmax_int))  
  
 print("outliers")
  # Evaluación de nodos internos 
  nodes <- outliers_int %>% pull(node) # extract all nodes of oulier intnodes 
  
  nodestocheck = data.frame(label=character(0), node=integer(0), ancestor=integer(0), 
                            edge.length=numeric(0), node.type=character(0), Outlier=integer(0), Filter=character(0)) 
  print("nodos extraidos")
  # We want to asses if some nodes are inside other nodes 
  # to remove common sites of VCF evaluation
  for (i in nodes){ # for each node
    ntips <- phangorn::Descendants(tr, i, "tips") # we obtain the nodes descendants tips
    
    temp <- (tr4df %>% filter(node %in% ntips[[1]]) %>% # filter the tree df to keep only descendants
               mutate(Outlier = !!i))  # add column with node name
    
    nodestocheck = rbind(nodestocheck, temp)
  }
  print("nodos comprovados")
  nodestocheck_def <- # check if any descendant tip of the node has to be tested
    nodestocheck %>%
    mutate(Filter = # Evaluate if there is a tip of node within the node that needs to be tested 
             case_when(( node.type == "internal" & edge.length > sdmax_int ) ~ "Test",
                       ( node.type == "tip" & edge.length > sdmax_tip ) ~ "Test",
                       TRUE ~ "Pass")) %>%
    full_join(., outliers_tips)  %>% # add outlier tips
    group_by(label) %>% 
    slice(which.max(Outlier)) %>% # remove duplicate tips if in two nodes
    group_by(Outlier) %>% # count tips in nodes
    add_count() 
  print("nodes def")
  if (nrow(nodestocheck_def) != 0){ # si hay alguna muestra que testar
    
    # Evaluo los aa consecutivos de cada tip de los nodos 
    
    outliersdef <-  nodestocheck_def %>% select(Outlier) %>% distinct() %>% pull() 
    print("outliersdef")
    for (outli in outliersdef) { #para cada nodo/tip
      teting_nodedf <- nodestocheck_def %>% filter( outli == Outlier)
      
      df5 <- vcf_list %>% 
        #uno la info del vcf con la de los nodos
        left_join(., teting_nodedf, by=c("isolate"="label")) %>%
        filter(REF != ALT) %>% # optimizacion
        #Eliminamos las posiciones comunes al clado
        group_by(Outlier) %>%
        distinct(POS, REF, ALT, .keep_all = T) %>%
        #Me quedo solo con los de testar
        filter(Filter == "Test") %>% 
        group_by(isolate) %>%
        mutate(GROUP = paste(isolate,cumsum(c(1L, diff(POS) > 1)), sep="-")) %>%
        group_by(isolate, GROUP) %>%
        #aQUI NECESITO TODA UNA EXPLICACION
        # diff(POS) TE LA LA DIFERENCIA ENTRE LA POS X Y X-1 (PORQUE EL LAG PREDETERMINADO ES 1, SE PUEDE CAmBIAR)
        # diff(POS) > 1 DEVUELVE TRUE/FALSE SI ENTRE UNA POSICION Y OTRA HAY MÁS DE 1
        # el 1L DE c(1L, diff(POS) > 1), CONVIERTE EL T/F A BINARIO
        # Y CUMSUM SUMA CUMULATIVAMENTE Y COMO ES BONARIO SE QUEDAN GRUPOS CONSECUTIVOS.
        filter(any(REF != "X") & any(ALT != "X") ) %>% # elimino todos los grupos que son gaps
        #group_by(GROUP, ALT) %>%
        filter(n() >= 3 ) %>%
        tally(ALT != "X") %>%
        filter(n >= 3 ) %>%
        group_by(isolate) %>%
        distinct(isolate) %>% 
        filter(isolate != "")
      
      removedf <- rbind(removedf, df5)
    }
  }
}

print("to names")
names <- anti_join(vcf_list %>% distinct(isolate), removedf) %>% 
  pull() %>%
  str_remove_all(., "^X") %>%
  str_replace_all(.,"P\\.", "P-")

write.table(names, args[3],  row.names = F,  col.names=F, quote = F)
write.table(removedf, args[4],  row.names = F,  col.names=F, quote = F)
