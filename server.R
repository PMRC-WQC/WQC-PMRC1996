#################################################################################################################
# Program for the PMRC 1996 G&Y model 
# Author: Cristian R. Montes
# Date  : January 2020
# Version 1.5
# Further distribution from the source code only with the author's permission.
#################################################################################################################

library(shiny)

#This is the definition for the set of equations out of Pienaars report

QMD  <- function(G, N)
{ #QMD stands for quadratic mean diameter
  #G   plot basal area
  #    N trees per acre  
  return(sqrt(G/(N*0.005454)))
}

#Product Definition
pulpwood  <- c(3,4)
sawtimber <- c(6,8)
veneer    <- c(8,10)


##################################################################################################################
##################################################################################################################
shinyServer(function(input, output, session) {

my_base <- reactive({

  if(input$Location == "Piedmont")            source("loblolly_Piedmont.R")
  if(input$Location == "Lower Coastal Plain") source("loblolly_LCP.R")
  if(input$Location == "Upper Coastal Plain") source("loblolly_UPC.R")
  idens = input$ThinDens #Density from the first thinning
  iAge  = input$Age[1]
  isi   = ifelse(input$SI != "Site Index",
                 Hi <- SI25(input$Hi, iAge),
                 Hi <- input$Hi)       #Site index from the unthinned counterpart

#Calculates everything for the unthinned scenario
  Ages.Base = seq(1, 35, 1)
  Ns.Base   = c(idens, N2(N1 = idens, 
                          SI = isi, 
                          A1 = 1, 
                          A2 = Ages.Base[-1]))
  Hs.Base   = H(isi, Ages.Base)
  PHWD      = 0
  
  Gs.Base   = G(H = Hs.Base, N = Ns.Base, A = Ages.Base, Nt = 0,  Nb = Ns.Base, At = 0, PHWD)
  
  Vs.Base   = V(H = Hs.Base, N = Ns.Base, G = Gs.Base,    A = Ages.Base)
  DWs.Base  = DW(H = Hs.Base, A = Ages.Base, N = Ns.Base, G = Gs.Base)
  Dq.Base   = QMD(G = Gs.Base, N = Ns.Base)

  Vpulp     = Vprod(V = Vs.Base, t = pulpwood[1],  D = Dq.Base,  N = Ns.Base, d = pulpwood[2])
  Vsaw      = Vprod(V = Vs.Base, t = sawtimber[1], D = Dq.Base,  N = Ns.Base, d = sawtimber[2])
  Vveneer   = Vprod(V = Vs.Base, t = veneer[1],    D = Dq.Base,  N = Ns.Base, d = veneer[2])
  
  Vpulp     = Vpulp - (Vsaw + Vveneer)
  Vsaw      = Vsaw  - Vveneer
  
  my_base   = data.frame(A       = Ages.Base,
                         N       = Ns.Base,
                         H       = Hs.Base,
                         G       = Gs.Base,
                         V       = Vs.Base,
                         DW      = DWs.Base,
                         Dq      = Dq.Base,
                         Vpulp   = Vpulp,
                         Vsaw    = Vsaw,
                         Vveneer = Vveneer,
                         Gpulp   = Vpulp/128 *  2.78,
                         Gsaw    = Vsaw/128 * 2.78,
                         Gveneer = Vveneer/128 * 2.78)
})

##Calculate Growth and Yield ###################################################
  calc_table<- reactive({
    
    if(input$Location == "Piedmont")           source("loblolly_Piedmont.R")
    if(input$Location =="Lower Coastal Plain") source("loblolly_LCP.R")
    if(input$Location =="Upper Coastal Plain") source("loblolly_UPC.R")
    
    iAge      = input$Age[1]
    eAge      = input$Age[2]
    thAge     = input$ThAge
    thDens    = input$ThinDens
    pulpwood  = input$specPr1
    sawtimber = input$specPr2
    veneer    = input$specPr3
    PHWD      = input$PHWD

    ifelse(input$SI != "Site Index",
          Hi <- input$Hi,
          Hi <- H(input$Hi, iAge))
    
    SI      = SI25(input$Hi, iAge)
    iDens   = input$N0 
    At      = thAge
    Nt      = thDens
    Nb      = N2(N1 = input$N0, SI = SI, A1 = iAge, A2 = thAge)
    
    Ages    = seq(iAge, eAge, 1)

        #Validate Consistency on the inputs 
    if (iAge>= eAge)
      eAge = iAge + 1
    
    if (thAge<= iAge)
      thAge = iAge + 1

    #Calculate Stand Values Before Thinning
    mx   =  round(N2(N1 = iDens, SI= SI, A1 = iAge, A2 = thAge),0)
    
    if(input$thin_type == "None")
    {
    Hs   = c(Hi,
             H2(H1 = Hi, A1 = iAge, A2 = Ages[-1]))
    Ns   = c(iDens, 
             N2(N1 = iDens, SI = SI, A1 = iAge, A2 = Ages[-1]))
    
    if(input$BA0 == 0)
      Gs   = G(H = Hs, N = Ns, A = Ages, Nt = 0, Nb = Ns, At = 0, PHWD)      
    else
    { Gs   = G2(G1 = input$BA0,
                A2 = Ages,
                A1 = iAge,
                H2 = Hs,
                H1 = Hs[1],
                N2 = Ns,
                N1 = Ns[1], 
                Nt = 0, 
                Nb = Ns,
                At = max(Ages)+1)
    }
    
    Vs   = V(H=Hs, N = Ns, G = Gs, A = Ages)
    DWs  = DW(H = Hs, A = Ages,N = Ns,G = Gs)
    
    As   = Ages
    }
    else
    {

    Ai    = seq(iAge, thAge, 1)
    Af    = seq(thAge, eAge, 1)  
    Ai.l  = length(Ai)
    Ai.L  = Ai.l - 1
    Af.l  = length(Af)
      
    Hs    = c(Hi,
              H2(H1 = Hi, A1 = iAge, A2 = Ai[-1]),
              H2(H1 = Hi, A1 = iAge, A2 = Af))
    Ns   = c(iDens,
             N2(N1 = iDens,  SI = SI, A1 = iAge,  A2 = Ai[-1]),
             N2(N1 = thDens, SI = SI, A1 = thAge, A2 = Af))
    Gi   = G(H = Hs[1],    N = Ns[1],      A = Ages[1],  Nt = 0, Nb = Ns[1],    At = 0, PHWD)
    Ga   = G(H = Hs[Ai.l], N = Ns[Ai.l+1], A = Ai[Ai.l], Nt = 0, Nb = Ns[Ai.l], At = 0, PHWD)
    Gs   = c(Gi,
             G2(G1 = Gi, 
                A2 = Ai[-1],
                A1 = iAge,
                H2 = Hs[2:Ai.l],
                H1 = Hs[1],
                N2 = Ns[2:Ai.l],
                N1 = Ns[1],
                Nb = Ns[2:Ai.l],
                Nt = Ns[2:Ai.l], 
                At = Ai),
             G(Hs[Ai.l],Ns[Ai.l+1],Ages[Ai.l],0,Ns[Ai.l],0, PHWD),
             G2(G1 = Ga,
                A2 = Af[-1],
                A1 = Ai[Ai.l],
                H2 = Hs[(Ai.l+2):(Af.l+Ai.l)],
                H1 = Hs[Ai.l+1],
                N2 = Ns[(Ai.l+2):(Af.l+Ai.l)],
                N1 = thDens,
                Nt = Ns[(Ai.l+1)] - thDens,
                Nb = Ns[Ai.l+1],
                At = thAge))
    
    Vs   = V(H = Hs, N = Ns, G = Gs, A = c(Ai, Af)) 
    DWs  = DW(H = Hs, A = c(Ai, Af),N = Ns,G = Gs)
    As  = c(Ai, Af)
    }
    
    Dq.s = QMD(G = Gs, N = Ns)
    
    Vpulps     = Vprod(V = Vs, t = pulpwood[1],  D = Dq.s, N = Ns, d = pulpwood[2])
    Vsaws      = Vprod(V = Vs, t = sawtimber[1], D = Dq.s, N = Ns, d = sawtimber[2])
    Vveneers   = Vprod(V = Vs, t = sawtimber[1], D = Dq.s, N = Ns, d = veneer[2])
  
    Vpulps     = Vpulps - (Vsaws + Vveneers)
    Vsaws      = Vsaws  - Vveneers
    
    Vpulps     = ifelse(Vpulps   >= 0, Vpulps,   0)
    Vsaws      = ifelse(Vsaws    >= 0, Vsaws,    0)
    Vveneers   = ifelse(Vveneers >= 0, Vveneers, 0)
    
    
    #Correct interface to reflect changes in density at the age of first thinning
    updateSliderInput(session, "ThAge",    min = input$Age[1], max = input$Age[2] - 1)
    updateSliderInput(session, "ThinDens", max = mx )
    updateSliderInput(session, "cstDistributionAge", min = input$Age[1], max = input$Age[2]-1)
   
    updateSliderInput(session, 
                      inputId = "Hi",
                      label = paste("Stand ", input$SI, " (ft)", sep =""),
                      value = input$Hi)
    
    
    my_table <- data.frame(A = As, N = Ns, H =  Hs, G = Gs, V = Vs, DW = DWs, Dq = Dq.s,
                           Vpulp   = Vpulps,
                           Vsaw    = Vsaws,
                           Vveneer = Vveneers,
                           Gpulp   = Vpulps/128 *  2.78,
                           Gsaw    = Vsaws/128 * 2.78,
                           Gveneer = Vveneers/128 * 2.78)
  
    
  })
  

##Volume extraction Calculation ################################################
  extraction <-reactive({
  
           my_table <- calc_table()         
           
           ln_inputs <- length(my_table$Vpulp)
           ln_initage<- my_table$A[1]
           ln_thage  <- input$ThAge
  
           
           Pulp = my_table$Vpulp[ln_inputs]
           Saw  = my_table$Vsaw[ln_inputs]
           Veneer = my_table$Vveneer[ln_inputs]
           if(input$chk1 == TRUE)
           {
            Pulprem    <- my_table$Vpulp[ln_thage - ln_initage +1] - my_table$Vpulp[2 + ln_thage - ln_initage] 
            Sawrem     <- my_table$Vsaw[ln_thage - ln_initage + 1] - my_table$Vsaw[2 + ln_thage - ln_initage]
            Veneerrem  <- my_table$Vveneer[ln_thage - ln_initage + 1] - my_table$Vveneer[2 + ln_thage - ln_initage]
           }
          else
          {
            Pulprem   = 0
            Sawrem    = 0
            Veneerrem = 0
          }  
          my_removals <- c(Pulprem = Pulprem, Sawrem = Sawrem, Veneerrem = Veneerrem, Pulp = Pulp, Saw = Saw, Veneer = Veneer)
          })
  

## Outputs #####################################################################
 output$SIPlot <- renderPlot({        my_base <- my_base()
                                      plot(c(0,35), c(0,100),
                                           type = "n",
                                           xlim = c(0,35),
                                           ylim = c(0,100),
                                           xlab = "Age (years)",
                                           ylab = "Dom. Height (ft)")
    
                                      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
                                           col = "lightblue")

                                      abline(v=seq(0,35,5), col = "white", lty = 2)
                                      abline(h=seq(0,100,10), col = "white", lty = 2)
                                    
        
                                      lines(H~A,
                                            data = my_base,
                                            lwd = 1,
                                            lty = 2,
                                            col ="darkgrey")
    
                                      new_data <- calc_table()
    
                                      lines(H ~A,
                                            data = new_data,
                                            lwd = 3,
                                            col = "red")
                                      
                                      SI = round(H2(new_data$H[1], A1 = new_data$A[1], 25),0)
                                      
                                      text(5,85.5,paste("Site Index ", SI), cex = 2.5, col = "white")
                                      text(5,86,paste("Site Index ", SI), cex = 2.5, col = "navy")
                                      
                                      with(new_data,
                                           points(c(A[1], A[length(A)]), 
                                                  c(H[1],H[length(H)]),
                                                  pch = 21,
                                                  col = "red",
                                                  bg = "yellow"))})

 output$SurvivalPlot <- renderPlot({  my_base <- my_base()
                                      plot(c(0,35),c(0,1100), type = "n",
                                           xlab = "Age (years)",
                                           ylab = "Survival (trees/acre)")
    
    
                                      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
                                          "lightblue")
                                          abline(v=seq(0,35,5), col = "white", lty = 2)
                                          abline(h=seq(0,1000,100), col = "white", lty = 2)
    
                                      lines(N ~ A,
                                            my_base,
                                            lty = 2,
                                            col ="darkgrey")
    
                                      new_data <- calc_table() 
   
                                      lines( N ~ A,
                                            new_data,
                                            lwd = 3,
                                            col = "red")
    
                                      with(new_data, points(c(A[1], A[length(A)]),
                                                    c(N[1],N[length(N)]),
                                                    pch = 21,
                                                    col = "red",
                                                    bg = "yellow"))})
  
 output$BasalArea <- renderPlot({my_base <- my_base()
                               plot(c(0,35),c(0,300), type = "n",
                               xlab = "Age (years)",
                               ylab = "Basal Area (ft2/acre)")
    
                                # Now set the plot region to grey
                               rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
                                   "lightblue")
                               abline(v=seq(0,35,5), col = "white", lty = 2)
                               abline(h=seq(0,500,50), col = "white", lty = 2)
    
                               lines(G ~ A, my_base, lty = 2,
                                     col ="darkgrey")
    
                               newData <- calc_table()
                               lines(G ~ A, newData, lwd = 3, col = "red")
                               with(newData, points(c(A[1], A[length(A)]), c(G[1], G[length(G)]), pch = 21, col = "red", bg = "yellow"))})
  
 output$Volume <- renderPlot({my_base <- my_base()
                              plot(c(0,35),c(0,10000), type = "n",
                                   xlab = "Age (years)",
                                   ylab = "Volume (ft3/acre)",
                                   main = "Volume ib")
    
                              rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
                                    "lightblue")
                              abline(v=seq(0,35,5), col = "white", lty = 2)
                              abline(h=seq(0,20000,2000), col = "white", lty = 2)
    
                              new_data <- calc_table()
    
                              lines(V ~ A, my_base, lty = 2,
                                    col = "darkgrey")
                              lines(V ~ A, new_data, col = "red", lwd = 3)
                              
                              lines(Vpulp ~ A , new_data, col = "blue",   lwd = 3)
                              lines(Vsaw ~ A, new_data,  col = "green",    lwd = 3)
                              lines(Vveneer ~ A, new_data, col = "brown", lwd = 3)
                              with(new_data, 
                                  points(c(A[1], A[length(A)]), c(V[1], V[length(V)]), pch = 21, col = "red", bg = "yellow"))
                              legend("topleft", legend=c("Base","Total", "Pulp", "Chip n Saw", "Sawtimber"),
                                     border ="navy",
                                     lty = c(2,1,1,1,1),
                                     lwd = c(1,3,3,3,3), col = c("darkgrey", "red", "blue", "green", "brown"))})
 
# output$VolGrowth<- renderPlotly({   my_base <- my_base()
#                                     new_data <-calc_table()
#                                     maxY <- 600#max(diff(new_data$V))
#                                     minY <- 0
#                                     p<- plot_ly(new_data, y= diff(V), x= A[-1],
#                                                 lines = list(shape ="spline"),  
#                                                 ylab  = "Volume growth (ft3/acre/year)",
#                                                 xlab  = "Age (years)",
#                                                 name = "Volume growth") 
#                                     layout(p, xaxis = list(title = "Age (years)",
#                                            range = c(0, 35),
#                                            autorange = F,
#                                            autotick = T),
#                                            yaxis = list(title = "Volume growth (ft3/acre/year)",
#                                                         range = c(minY, maxY)))
#                                  })
 
 output$DryWeight <- renderPlot({my_base <- my_base()
                                 plot(c(0,35),c(0,100), type = "n",
                                      xlab = "Age (years)",
                                      ylab = "Dry Weight (Ton/acre)",
                                      main = "Total Weight ib.")
                                   
                                 # Now set the plot region to grey
                                 rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
                                        "lightblue")
                                 abline(v=seq(0,35,5), col = "white", lty = 2)
                                 abline(h=seq(0,100,10), col = "white", lty = 2)
                                   
                                 new_data <- calc_table()
                                   
                                 lines(DW ~ A, my_base, lty = 2,
                                       col = "darkgrey")
                                 lines(DW ~ A, new_data, col = "red", lwd = 3)
                                 with(new_data, 
                                      points(c(A[1], A[length(A)]), c(DW[1], DW[length(DW)]), pch = 21, col = "red", bg = "yellow"))})  
  
 output$QMD         <- renderPlot({ my_base <- my_base()
                                  plot(c(0,35),c(0,15), type = "n",
                                         xlab = "Age (years)",
                                         ylab = "Quadratic Mean Diameter (in)")
   
                                  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
                                        "lightblue")
                                        abline(v=seq(0,35,5), col = "white", lty = 2)
                                        abline(h=seq(0,25,5), col = "white", lty = 2)
   
                                        new_data <- calc_table()
   
                                        lines(Dq ~ A, my_base, lty = 2,
                                              col = "darkgrey")
                                        lines(Dq ~ A, new_data, col = "red", lwd = 3)
                                  with(new_data, 
                                        points(c(A[1], A[length(A)]), 
                                               c(Dq[1], Dq[length(DW)]),
                                               pch = 21, col = "red", bg = "yellow"))})
 
 output$Extraction   <- renderPlot({my_base <- my_base()
                             my_cols <- c("Gold", "darkseagreen1", "plum2")
                             my_table <- calc_table()
                             outs <- matrix(extraction(),3,2)
                             
                             maxVol <- max(my_table$V[length(my_table[,1])], sum(outs[,1]))
                             
                             colnames(outs) <- c("Thinning Removals","Harvest Removals")  
                             rownames(outs) <- c("Pulp", "Chip n Saw", "Sawtimber")
                             barplot(outs, col = my_cols,
                                     ylab = "Volume (ft3/acre)",
                                     ylim = c(0,maxVol),
                                     border = "grey")
                             rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
                                    "lightblue")
                             
                             abline(h = seq(0,maxVol, 500), lty = 2, col = "white")
                             
                             barplot(outs, col = my_cols,
                                     ylab = "Volume (ft3/acre)",
                                     ylim = c(0,maxVol),
                                     border = "grey",
                                     inside = FALSE,
                                     add = TRUE)

                             box(lwd = 3)
                             legend("topleft", 
                                    fill = my_cols,
                                    legend = c("Pulp", "Chip n Saw", "Sawtimber"),
                                    border = "grey",
                                    density = 100)
                             
                             })
 
 output$Yield.Table <- renderTable({
                                   Yield.table <- calc_table()                               
                                   colnames(Yield.table)<-c("Age (years)", 
                                                            "Dens. (trees/acre)", 
                                                            "Height (ft)", 
                                                            "BArea (ft2/acre)", 
                                                            "Vol ob (ft3/acre)",
                                                            "Dry Weight (ton/acre)",
                                                            "Dq (inches)",
                                                            "Pulp (ft3/acre)",
                                                            "Chip-n-Saw (ft3/acre)",
                                                            "Sawtimber (ft3/acre)",
                                                            "GW Pulp (ton/acre)",
                                                            "GW Chip-n-saw (ton/acre)",
                                                            "GW Sawtimber (ton/acre)")
                                   Yield.table
                                   },include.rownames=FALSE)
 

#output$Marginal  <- renderPlotly({
#                                financial <-calc_costs()
#                                maxY <- max(max(diff(financial$NPV)),max(diff(financial$BLV)))
#                                minY <- min(min(diff(financial$NPV), min(diff(financial$BLV))))
#                                p<- plot_ly(financial, y= diff(NPV), x= Age[-1],
#                                     lines = list(shape ="spline"),  
#                                     ylab = "Marginal Gain ($/acre)",
#                                     xlab = "Age (years)",
#                                     ylim = c(minY, maxY),
#                                     name = "delta NPV") %>%
#                                add_trace(y = diff(BLV), name = "delta BLV", line = list(shape = "spline"))
                                
#                                layout(p, xaxis = list(title = "Age (years)",
#                                                      range = c(min(financial$Age), max(financial$Age)),
#                                                      autorange = F,
#                                                     autotick = T),
#                                          yaxis = list(title = "Marginal Gain ($/acre)"))
#})
  

output$Distribution <- renderPlot({
                          Yield.table <- calc_table()  
                          Age         <- input$cstDistributionAge
                          Yield.table <- Yield.table[Yield.table$A == Age,]
                          S           <- H2(Yield.table$H, 25, Age)
                          values      <- GetDiameterDistribution(Age = Age,
                                                                 N   = Yield.table$N,
                                                                 H   = Yield.table$H,
                                                                 G   = Yield.table$G,
                                                                 S   = S)
                          barplot(values, names = seq(1:20), ylim = c(0,500), ylab = "Class population (trees/acre)", xlab = "Diameter Class (in)")
                          box()
})

output$downloadData <- downloadHandler(
  filename = function() {
    paste(input$dataset, ".csv", sep = "")
  },
  content = function(file) {
    write.csv(calc_table(), file, row.names = FALSE)
  }
)
})