
server <- shinyServer(function(input, output, session) {

  dat <- reactiveVal()
  groups <- reactiveVal(data.frame(matrix(nrow=0,ncol=3)))
  btwn_groups <- reactiveVal(data.frame(matrix(nrow=0,ncol=4)))
  contrast_res <- reactiveVal()
  tmp_levels <- c("SummerActive","30C","25C","20C","12C","4C")
  target<- reactiveVal()
  iformula <- reactiveVal()
  ig <- 1
  ibg <- 1
  
  
  ###
  #def within_string
  within_string <- function(igroup){
    membs <- groups %>% 
      filter(Groups==igroup) %>% 
      pull(Members)
    rel <- groups %>% 
      filter(Groups==igroup) %>% 
      pull(WithinRel) %>% unique()
    print(rel)
    if (rel=='Difference'){
      istr <- paste0('(',paste(paste0(rep(input$itiss,length(membs)),'.',membs), collapse =  ' - '),')')

    } else if (rel=='Average'){
      istr <- paste0('(', paste(paste0(rep(input$itiss,length(membs)),'.',membs), collapse =  ' + '), ')','/',length(membs))
    }
    return(istr)
  }
  
  #####
  #load rds
  observeEvent(input$load_results, {
    req(input$datfile)
    
    dat <<- readRDS(input$datfile$datapath)
    targets <<- read_csv('../data/deseq/all_tiss_meta.csv') %>% 
      mutate(group=paste(tissue,temp,sep='.')) %>% 
      filter(sample %in% colnames(dat$y)) %>% 
      mutate(files=sample) %>% 
      mutate(color=case_when(tissue=='kidney'~'orange',
                             tissue=='brain' ~ 'purple',
                             tissue=='heart' ~ 'red', 
                             tissue=='liver' ~ 'blue'))
    
    })

  #####
  #add group
  observeEvent(input$add_group,{
    membs <- input$group_in
    
    if (ig<2){groups<<- data.frame('Groups'=paste('G',ig,sep=''), 
                                             'Members'=membs,
                                             'WithinRel'=input$within_rel) 
    } else { groups <<- rbind(groups, 
               data.frame('Groups'=paste('G',ig,sep=''), 
                          'Members'=membs,
                          'WithinRel'=input$within_rel))
    }
    
    ig <<- ig + 1
    
    #update group table
    output$groups<- DT::renderDataTable({
      DT::datatable(groups, list(pageLength = 10, scrollX=T))
    })
    
    #update group lists
    output$group_select <- renderUI({
      req(ig>1)
      selectInput('gedit_select','Select Group',
                  choices = unique(groups$Groups), selected = 'G1')
    })
    
    output$group2_select <- renderUI({
      req(ig>2)
      selectInput('gedit2_select','Select Group 2',
                  choices = unique(groups$Groups), selected='G2')
    })
  })
  
  

  #####
  # generate final formula, plot graph

  observeEvent(input$run_contr, {
    #update table
    req(input$gedit_select)
    req(input$gedit2_select)
    istr1<-within_string(input$gedit_select)
    istr2<-within_string(input$gedit2_select)
    btwn_groups <<- data.frame('Group1'=istr1,
                               'Group2'=istr2,
                               'Formula'=paste0(istr1,' - ', istr2),
                               'Sign'=input$tsign)
    
    iformula <<- btwn_groups %>% pull(Formula) %>% unique()
    
    #output$formula <- renderText({iformula})
    
    
    #update graph
    output$model_graph <- renderPlot({
      
      print(iformula)
      contrast.matrix <- makeContrasts(contrasts=iformula, levels= dat$design )
      print(typeof(input$ngenes))
      ifit <- contrasts.fit(dat$fit, contrast.matrix)
      ifit <- eBayes(ifit)
      iresults <- decideTests(ifit)
      if (input$tsign=='+'){
        out <- topTable(ifit, p.value=0.01, number = as.integer(input$ngenes) ) %>% filter(t>0) %>% rownames_to_column('gene')
        subtit <- 't > 0'
      } else {
        out <- topTable(ifit, p.value=0.01, number = as.integer(input$ngenes) ) %>% filter(t<0) %>% rownames_to_column('gene')
        subtit <- 't < 0'
      }
      contrast_res <<- topTable(ifit, p.value=0.01, number = as.integer(input$ngenes)*2) %>% rownames_to_column('gene')
      glist <- out %>% pull(gene)
      lsamps <- targets %>% filter(tissue==input$itiss)
      icolor <- lsamps %>% pull(color) %>% unique()
      esub <- dat$y[glist,]$E
      exprdat <- data.frame(esub) %>% 
        rownames_to_column('gene') %>% 
        pivot_longer(cols=-c(gene), names_to='sample', values_to = 'expr')
      gdat <- merge(lsamps, exprdat, by='sample') %>% 
        mutate(temp=factor(temp, levels=c('4C','12C','20C','25C','30C','SummerActive'))) %>% 
        group_by(gene, temp) %>% summarise(avg_expr=mean(expr))
      ggplot(gdat,aes(x=temp, y=avg_expr, group=gene, label=gene))+
        geom_line(color=icolor)+
        #geom_label_repel(aes(label=gene), nudge_x = 1, label.size = 0.05)+
        theme_bw()+
        ggtitle(iformula,subtitle = subtit)
      
      
    })
    
    #update table
    output$contrast_results <- DT::renderDataTable({
      out <- contrast_res 
      if (input$tsign=='+'){
        out <- out %>% filter(t>0)
      } else {
        out <- out %>% filter(t<0)
      }
      DT::datatable(out, list(pageLength = 10, scrollX=T))
    },
    options = list(
      autoWidth = TRUE,
      columnDefs = list(list(width = '100px', targets = "_all"))
    ))
  })
  
  
  
  #####
  #render model plot
  
  
    
  
    
  #####
  #Refresh table
    

  
})