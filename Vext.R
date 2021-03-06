#------------------------------------------C�LCULO DE VREF POR DIVERSOS M�TODOS-------------------------------------------------------------------

#Par�metros:
#velocidad: vector velocidad de d�nde se extraer� la serie de m�ximos
#fecha: vector fecha en formato "AAAA-MM-DD" (cualquier separador v�lido) para Periodic Maxima.
#TR: Per�odo de retorno (en a�os). 50 por defecto.
#method: M�todo para ajustar Weibull para EWTS's. Puede ser: "Mom", "LS", "Wasp", "MLE" � "MMLE".
#equation: M�todo para obtener Vref. Puede ser: "EWTS_E", "EWTS_G", "EWTS_D", "Tormentas","Periodic_M�xima" � "Gumbel"
#dias: separaci�n de d�as para localizar tormentas independientes (con funci�n MIS)
#tipo_per: Para la ecuaci�n "Periodic_Maxima". Puede ser "A�o" o "Mes" y buscar� los m�ximos entre estos per�odos.
#l�mite: m�nimo para el m�todo de las tormentas. Se eliminan los m�ximos inferiores a este valor.
#prob_method: Para Tormentas y Periodic_Maxima. Dos formas distintas para calcular la probabilidad a partir del n�mero de tormentas.
 #Usar 1 o 2.
#dias_sep: separaci�n en d�as entre los m�ximos detectados. Tiene que ser menor que "dias".
#gumbel_method: Puede ser "moments" � "MLE". Para estimar los par�metros por uno de esos dos m�todos.
#diezmin: TRUE o FALSE. Para el m�todo de las tormentas y periodic maxima. TRUE para 10minutales, FALSE para segundales

Vext<-function(velocidad, fecha=NULL,Tr=50, method="Mom", equation, dias=4,
               tipo_per="A�o",limite=19, prob_method=1, dias_sep=2, gumbel_method, diezmin=TRUE){
  
  n=23037
  
  if(is.vector(velocidad)==FALSE){
    stop("velocidad debe ser un vector")
  }
  
  if(length(which(velocidad>90))!=0){
    stop("velocidad debe ser filtrado. Valores superiores a 90 encontrados")
  }
  
  if(is.numeric(velocidad)==FALSE){
    stop("velocidad debe ser num�rico")
  }
  
  if(is.character(method)==FALSE){
    stop("method es un argumento inv�lido. Debe ser 'character'")
  }
  
  if(is.character(gumbel_method)==FALSE){
    stop("gumbel_method es un argumento inv�lido. Debe ser 'character'")
  }
  
  if(gumbel_method!="moments" & gumbel_method!="MLE"){
    stop("gumbel_method debe ser 'moments' � 'MLE'")
  }
  
  if(method!="LS" & method!="Mom" & method!="Wasp" & method!="MLE" & method!="MMLE"){
    stop("method es un argumento inv�lido. Debe ser: 'LS', 'Mom', 'MLE', 'MMLE' � 'Wasp'")
  }
  
  if(is.numeric(dias)==FALSE){
    stop("dias debe ser num�rico")
  }
  
  if(is.numeric(dias_sep)==FALSE){
    stop("dias_sep debe ser num�rico")
  }
  
  if(dias<dias_sep){
    stop("'dias' debe ser mayor que 'dias_sep'")
  }
  
  if(is.numeric(limite)==FALSE){
    stop("limite debe ser num�rico")
  }
  
  if(is.numeric(prob_method)==FALSE){
    stop("prob_method debe ser num�rico. Puede ser 1 � 2.")
  }
  
  if(tipo_per!="A�o" & tipo_per!="Mes"){
    stop("tipo_per tine que ser 'A�o' o 'Mes' ")
  }
  
  if(is.numeric(Tr)==FALSE){
    stop("Tr tiene que ser num�rico")
  }
  
  WeibullParam<-function(V1,k=seq(0.1,10,0.01), method){
    
    library(ggthemes)
    library(extrafont)
    library(ggplot2)
    library(gridExtra)
    
    barfill <- "#4271AE"
    barlines <- "#1F3552"
    
    V3<-floor(V1)
    
    if(is.vector(V1)==FALSE){
      stop("V1 debe ser un vector")
    }
    
    if(length(which(V1>90))!=0){
      stop("V1 debe ser filtrado. Valores superiores a 90 encontrados")
    }
    
    if(is.numeric(V1)==FALSE){
      stop("V1 debe ser num�rico")
    }
    
    if(is.character(method)==FALSE){
      stop("method es un argumento inv�lido")
    }
    
    if(is.vector(k)==FALSE | is.numeric(k)==FALSE){
      stop("k debe ser un vector num�rico")
    }
    
    if(method=="LS"){
      
      ordenar<-order(V1)
      
      V1<-V1[ordenar]
      frecuencia<-seq(1,length(V1),1)/(length(V1)+1)
      
      y<-log(-log(1-frecuencia))
      x<-log(V1)
      
      rx<-lm(y~x)
      
      df1<-as.data.frame(matrix(nrow = length(x), ncol = 2)); df1$x<-x; df1$y<-y
      
      p6<-ggplot(data = df1, aes(x=x,y=y))+
        geom_point(size=3, colour="red", fill=I("red"))+
        geom_smooth(method = lm, formula = y~x, lwd=1.5)+
        scale_x_continuous(name = "ln(U)",
                           limits=c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))) +
        scale_y_continuous(name = "ln(ln(1/1-F(u)))", limits = c(min(y, na.rm = TRUE),max(y, na.rm = TRUE))) +
        labs(colour="Ajuste")+
        ggtitle("Estimaci�n de par�metros por m�todo gr�fico") +
        theme_economist() +
        theme(legend.position = "bottom", legend.direction = "horizontal",
              legend.box = "horizontal",
              legend.key.size = unit(1, "cm"),
              plot.title = element_text(family="Bauhaus 93", size = 18),
              text = element_text(family = "Bauhaus 93"),
              axis.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 9),
              legend.title=element_text(face = "bold", size = 9, family = "Bauhaus 93"))
      
      k<-as.vector(rx$coefficients[2])
      intercept<-as.vector(rx$coefficients[1])
      A<-exp(-intercept/k)
      
      
      #C�lculo de Chi_Cuadrado y RMSE:
      
      V1<-V1[is.na(V1)==FALSE]
      
      frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
      
      for(i in 1:length(frecuencia)){
        frecuencia1[i]<-(k/A)*((V1[i]/A)^(k-1))*exp(-(V1[i]/A)^k)
      }
      
      suma_chi<-0
      
      for(i in 1:length(frecuencia1)){
        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
      }
      
      media<-mean(frecuencia, na.rm = TRUE)
      suma11<-0
      
      for(i in 1:length(frecuencia1)){
        suma11<-suma11+(frecuencia[i]-media)^2
      }
      
      R<-1-(suma_chi/suma11)
      
      Chi_Cuadrado<-suma_chi/(length(V1)-2)
      
      RMSE<-(suma_chi/length(V1))^0.5
      

      Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE,R_2=R, graph=p6)
      return(Coeficientes)
    }
    
    if(method=="Mom"){
      
      vmedia<-mean(V1[is.na(V1)==FALSE])
      smedia<-sd(V1[is.na(V1)==FALSE])
      
      numerador<-vector(mode = "numeric", length = length(k))
      denominador<-vector(mode="numeric", length = length(k))
      formula<-vector(mode = "numeric", length = length(k))
      
      for(j in 1:length(k)){
        numerador[j]<-(smedia^2)*(gamma(1+(1/k[j]))^2)
        denominador[j]<-(vmedia^2)*(gamma(1+(2/k[j]))-gamma(1+(1/k[j]))^2)
        formula[j]<-(numerador[j]/denominador[j])-1
      }
      
      minimo<-min(abs(formula))
      loc_min<-which(abs(formula)==minimo)
      
      k<-k[loc_min]
      A<-vmedia/gamma(1+(1/k))
      
      #C�lculo de Chi_Cuadrado y RMSE:
      
      Vmedia1<-mean(V1[is.na(V1)==FALSE])
      V4<-as.vector(table(V3))
      suma5<-sum(V4)
      frecuencia<-V4/suma5
      Vbinmedio<-as.vector(tapply(V1,V3,mean))
      
      frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
      
      for(i in 1:length(frecuencia)){
        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
      }
      
      suma_chi<-0
      
      for(i in 1:length(frecuencia1)){
        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
      }
      
      media<-mean(frecuencia, na.rm = TRUE)
      suma11<-0
      
      for(i in 1:length(frecuencia1)){
        suma11<-suma11+(frecuencia[i]-media)^2
      }
      
      R<-1-(suma_chi/suma11)
      
      Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)
      
      RMSE<-(suma_chi/length(Vbinmedio))^0.5
      
      
      
      z<-as.numeric(names(tapply(V1,V3,mean)))
      
      frecuencia1_df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3)); 
      
      colnames(frecuencia1_df)<-c("frecuencia", "z", "frecuencia1")
      
      frecuencia1_df$frecuencia<-frecuencia; frecuencia1_df$z<-z; 
      
      frecuencia1_df$frecuencia1<-frecuencia1
      
      
      
      par(mfrow=c(1,1))
      
      barfill <- "#4271AE"
      barlines <- "#1F3552"
      
      p7 <- ggplot(frecuencia1_df) +
        geom_col(aes(x =z+0.5, y=frecuencia),
                 colour = barlines, fill = barfill) +
        scale_x_continuous(name = "Velocidad (m/s)",
                           breaks = seq(1,length(frecuencia1),1),
                           limits=c(0, length(frecuencia))) +
        scale_y_continuous(name = "Frecuencia", limits = c(0,0.15)) +
        geom_line(data = frecuencia1_df, aes(x=z+0.5, y=frecuencia1, colour="red"), size=1.5)+
        labs(colour="Ajuste")+
        ggtitle("Distribuci�n ajustada a Weibull por Momentos") +
        theme_economist() +
        theme(legend.position = c(0.75,0.75), legend.direction = "horizontal",
              legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
              legend.box = "horizontal",
              legend.key.size = unit(1, "cm"),
              plot.title = element_text(family="comic-sans", size = 22),
              text = element_text(family = "comic-sans"),
              axis.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 9),
              legend.title=element_blank())+
        scale_colour_discrete("Point",labels=c("Ajuste Weibull"))
      
      p7
      
      Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p7)
      

      return(Coeficientes)
    }
    
    if(method=="Wasp"){
      
      V1<-V1[is.na(V1)==FALSE]
      
      ordenar<-order(V1); V1<-V1[ordenar]
      
      V2<-abs(V1-mean(V1, na.rm = TRUE))
      
      loc<-which(V2==min(V2, na.rm = TRUE))
      
      X<-1-(loc[length(loc)]/(length(V1)+1))
      
      Vmedia<-mean(V1, na.rm = TRUE); suma<-sum(V1^3, na.rm = TRUE); N<-length(V1)
      
      
      #Aplicamos la f�rmula
      
      numerador1<-vector(mode = "numeric", length = length(k))
      numerador2<-vector(mode = "numeric", length = length(k))
      formula<-vector(mode = "numeric", length = length(k))
      
      for(j in 1:length(k)){
        numerador1[j]<-(-log(X))^(1/k[j])
        numerador2[j]<-(suma/(N*gamma(1+(3/k[j]))))^(1/3)
        
        formula[j]<-((numerador1[j]*numerador2[j])/Vmedia)-1
      }
      
      #Buscamos la k para la que el resultado est� m�s pr�ximo a 0:
      
      minimo<-min(abs(formula))
      loc_min<-which(abs(formula)==minimo)
      k<-k[loc_min]
      
      
      A<-(suma/(N*gamma(1+(3/k))))^(1/3)
      
      #C�lculo de Chi_Cuadrado y RMSE:
      
      Vbinmedio<-as.vector(tapply(V1, V3, mean))
      frecuencia1<-vector(mode = "numeric", length = length(Vbinmedio))

      for(i in 1:length(Vbinmedio)){
        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
      }

      suma_chi<-0

      for(i in 1:length(frecuencia1)){
        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
      }

      media<-mean(frecuencia, na.rm = TRUE)
      suma11<-0

      for(i in 1:length(frecuencia1)){
        suma11<-suma11+(frecuencia[i]-media)^2
      }

      R<-1-(suma_chi/suma11)

      Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)

      RMSE<-(suma_chi/length(Vbinmedio))^0.5

      z<-as.numeric(names(tapply(V1,V3,mean)))

      par(mfrow=c(1,1))

      frecuencia1_df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3)); colnames(frecuencia1_df)<-c("frecuencia", "z", "frecuencia1")

      frecuencia1_df$frecuencia<-frecuencia; frecuencia1_df$z<-z; frecuencia1_df$frecuencia1<-frecuencia1

      barfill <- "#4271AE"
      barlines <- "#1F3552"

      p7 <- ggplot(frecuencia1_df) +
        geom_col(aes(x =z+0.5, y=frecuencia),
                 colour = barlines, fill = barfill) +
        scale_x_continuous(name = "Velocidad (m/s)",
                           breaks = seq(1,length(frecuencia1),1),
                           limits=c(0, length(frecuencia))) +
        scale_y_continuous(name = "Frecuencia", limits = c(0,0.15)) +
        geom_line(data = frecuencia1_df, aes(x=z+0.5, y=frecuencia1, colour="red"), size=1.5)+
        labs(colour="Ajuste")+
        ggtitle("Distribuci�n ajustada a Weibull por Wasp") +
        theme_economist() +
        theme(legend.position = c(0.75, 0.75), legend.direction = "horizontal",
              legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
              legend.box = "horizontal",
              legend.key.size = unit(1, "cm"),
              plot.title = element_text(family="comic-sans", size = 22),
              text = element_text(family = "comic-sans"),
              axis.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 9),
              legend.title=element_blank())+
        scale_colour_discrete("Point",labels=c("Ajuste Weibull"))

      p7
      
      Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p7)

      return(Coeficientes)
    }
    
    if(method=="MLE"){
      
      N<-length(V1[is.na(V1)==FALSE])
      
      suma<-0
      
      for(i in 1:N){
        if(is.na(V1[i])==FALSE){
          suma<-suma+log(V1[i])
        }
      }
      
      suma1<-0
      
      for(i in 1:N){
        if(is.na(V1[i])==FALSE){
          suma1<-suma1+(V1[i])^k
        }
      }
      
      suma2<-0
      
      for(i in 1:N){
        if(is.na(V1[i])==FALSE){
          suma2<-suma2+(((V1[i])^k)*log(V1[i]))
        }
      }
      
      #Calculamos la funci�n:
      numerador1<-k*N*suma2
      numerador2<-k*suma1*suma
      denominador<-N*suma1
      formula<-((numerador1-numerador2)/denominador)-1
      
      #Obtenemos k en base al valor m�s pr�ximo a 0 obtenido en la f�rmula:
      minimo<-min(abs(formula))
      loc_min<-which(abs(formula)==minimo)
      k<-k[loc_min]
      
      
      #Obtenemos A a partir de k:
      suma4<-0
      
      for(i in 1:N){
        if(is.na(V1[i])==FALSE){
          suma4<-suma4+((V1[i])^k)
        }
      }
      
      A<-(suma4/N)^(1/k)
      
      #C�lculo de Chi_Cuadrado y RMSE:
      
      Vmedia<-mean(V1[is.na(V1)==FALSE])
      V4<-as.vector(table(V3))
      suma5<-sum(V4)
      frecuencia<-V4/suma5
      Vbinmedio<-as.vector(tapply(V1,V3,mean))
      
      frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
      
      for(i in 1:length(frecuencia)){
        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
      }
      
      suma_chi<-0
      
      for(i in 1:length(frecuencia1)){
        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
      }
      
      media<-mean(frecuencia, na.rm = TRUE)
      suma11<-0
      
      for(i in 1:length(frecuencia1)){
        suma11<-suma11+(frecuencia[i]-media)^2
      }
      
      R<-1-(suma_chi/suma11)
      
      Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)
      
      RMSE<-(suma_chi/length(Vbinmedio))^0.5
      
      z<-as.numeric(names(tapply(V1,V3,mean)))
      
      par(mfrow=c(1,1))
      
      frecuencia1_df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3)); colnames(frecuencia1_df)<-c("frecuencia", "z", "frecuencia1")
      
      frecuencia1_df$frecuencia<-frecuencia; frecuencia1_df$z<-z; frecuencia1_df$frecuencia1<-frecuencia1
      
      p7 <- ggplot(frecuencia1_df) +
        geom_col(aes(x =z+0.5, y=frecuencia),
                 colour = barlines, fill = barfill) +
        scale_x_continuous(name = "Velocidad (m/s)",
                           breaks = seq(1,length(frecuencia1),1),
                           limits=c(0, length(frecuencia))) +
        scale_y_continuous(name = "Frecuencia", limits = c(0,0.15)) +
        geom_line(data = frecuencia1_df, aes(x=z+0.5, y=frecuencia1, colour="red"), size=1.5)+
        labs(colour="Ajuste")+
        ggtitle("Distribuci�n ajustada a Weibull por M�xima Verosimilitud") +
        theme_economist() +
        theme(legend.position = c(0.75, 0.75), legend.direction = "horizontal",
              legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
              legend.box = "horizontal",
              legend.key.size = unit(1, "cm"),
              plot.title = element_text(family="Bauhaus 93", size = 18),
              text = element_text(family = "comic-sans"),
              axis.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 9),
              legend.title=element_blank())+
        scale_colour_discrete("Point",labels=c("Ajuste Weibull"))
      
      p7
      
      Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p7)
      
      
      return(Coeficientes)
    }
    
    if(method=="MMLE"){
      
      N<-length(V1[is.na(V1)==FALSE])
      
      suma<-0
      
      for(i in 1:N) {
        if(is.na(V1[i])==FALSE){
          suma<-suma+(log(V1[i])*log(V1[i]))
          
        }
      }
      
      suma1<-0
      
      for(i in 1:N){
        if(is.na(V1[i])==FALSE){
          suma1<-suma1+log(V1[i])
        }
      }
      
      numerador<-N*(N-1)
      denominador<-(N*suma)-((suma1)^2)
      k<-(pi/sqrt(6))*(numerador/denominador)^0.5
      
      suma3<-0
      
      for(i in 1:N){
        if(is.na(V1[i])==FALSE){
          suma3<-suma3+(V1[i])^k
        }
      }
      
      A<-(suma3/N)^(1/k)
      
      #C�lculo de Chi_Cuadrado y RMSE:
      
      Vmedia1<-mean(V1[is.na(V1)==FALSE])
      V4<-as.vector(table(V3))
      suma5<-sum(V4)
      frecuencia<-V4/suma5
      Vbinmedio<-as.vector(tapply(V1,V3,mean))
      
      frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
      
      for(i in 1:length(frecuencia)){
        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
      }
      
      suma_chi<-0
      
      for(i in 1:length(frecuencia1)){
        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
      }
      
      media<-mean(frecuencia, na.rm = TRUE)
      suma11<-0
      
      for(i in 1:length(frecuencia1)){
        suma11<-suma11+(frecuencia[i]-media)^2
      }
      
      R<-1-(suma_chi/suma11)
      
      Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)
      
      RMSE<-(suma_chi/length(Vbinmedio))^0.5
      
      z<-as.numeric(names(tapply(V1,V3,mean)))
      
      par(mfrow=c(1,1))
      
      frecuencia1_df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3)); colnames(frecuencia1_df)<-c("frecuencia", "z", "frecuencia1")
      
      frecuencia1_df$frecuencia<-frecuencia; frecuencia1_df$z<-z; frecuencia1_df$frecuencia1<-frecuencia1
      
      p7 <- ggplot(frecuencia1_df) +
        geom_col(aes(x =z+0.5, y=frecuencia),
                 colour = barlines, fill = barfill) +
        scale_x_continuous(name = "Velocidad (m/s)",
                           breaks = seq(1,length(frecuencia1),1),
                           limits=c(0, length(frecuencia))) +
        scale_y_continuous(name = "Frecuencia", limits = c(0,0.15)) +
        geom_line(data = frecuencia1_df, aes(x=z+0.5, y=frecuencia1, colour="red"), size=1.5)+
        labs(colour="Ajuste")+
        ggtitle("Distribuci�n ajustada a Weibull por Modificaci�n M�xima Verosimilitud") +
        theme_economist() +
        theme(legend.position = c(0.75, 0.75), legend.direction = "horizontal",
              legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
              legend.box = "horizontal",
              legend.key.size = unit(1, "cm"),
              plot.title = element_text(family="comic-sans", size = 16),
              text = element_text(family = "comic-sans"),
              axis.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 9),
              legend.title=element_blank())+
        scale_colour_discrete("Point",labels=c("Ajuste Weibull"))
      
      p7

      Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p7)
      
      
      return(Coeficientes)
    }
    
    else{
      stop("method es un argumento inv�lido")
    }
  }
  
  Maximos_MIS<-function(velocidad, dias, dias_sep, diezmin=TRUE){
    
    if(is.vector(velocidad)==FALSE){
      stop("velocidad tiene que ser un vector")
    }
    
    if(is.numeric(velocidad)==FALSE | is.numeric(dias)==FALSE){
      stop("Velocidad y dias tienen que ser num�ricos")
    }
    
    if(sum(which((velocidad>60)!=0))){
      stop("Velocidades mayores que 60 encontradas. Debe filtrarse el vector")
    }
    
    if(diezmin==TRUE){
      
      paso<-6*24*dias
      
      long<-floor(length(velocidad)/paso)
      
      velocidad<-velocidad[1:(long*paso)]
      
      x<-seq(1,length(velocidad),paso)
      
      Maximos_MIS<-vector(mode = "numeric", length = length(x))
      localizador<-vector()
      
      for(i in 2:length(x)){
        Maximos_MIS[i-1]<-max(velocidad[x[i-1]:x[i]], na.rm = TRUE)
        localizador[i-1]<-which(velocidad[x[i-1]:x[i]]==max(velocidad[x[i-1]:x[i]], na.rm = TRUE))
      }
      
      localizador<-localizador+(x[1:length(localizador)]-1)
      
      diff<-vector()
      
      for(i in 2:length(Maximos_MIS)){
        diff[i-1]<-localizador[i]-localizador[i-1]
      }
      
      loc<-which(diff>(6*24*dias_sep)); maximos<-vector(); maximos<-Maximos_MIS[loc]
      
      resultado<-list(Maximos=maximos, localizador=localizador)
      
      return(maximos)
    }
    
    else{
      paso<-3600*24*dias
      
      long<-floor(length(velocidad)/paso)
      
      velocidad<-velocidad[1:(paso*long)]
      
      x<-seq(1,length(velocidad),paso)
      
      Maximos_MIS<-vector(mode = "numeric", length = length(x))
      localizador<-vector()
      
      for(i in 2:length(x)){
        Maximos_MIS[i-1]<-max(velocidad[x[i-1]:x[i]], na.rm = TRUE)
        localizador[i-1]<-which(velocidad[x[i-1]:x[i]]==max(velocidad[x[i-1]:x[i]], na.rm = TRUE))
      }
      
      localizador<-localizador+(x[1:length(localizador)]-1)
      
      diff<-vector()
      
      for(i in 2:length(Maximos_MIS)){
        diff[i-1]<-localizador[i]-localizador[i-1]
      }
      
      loc<-which(diff>(3600*24*dias_sep)); maximos<-vector(); maximos<-Maximos_MIS[loc]
      
      resultado<-list(Maximos=maximos, localizador=localizador)
      
      return(maximos)
    }
    
  }
  
  Weibull<-WeibullParam(velocidad, method = method)
  
  k<-Weibull$k
  
  vmedia<-mean(velocidad, na.rm = TRUE)
  
  if(equation=="EWTS_E"){
    
    a<-1/gamma(1+(1/k))
    
    b<-exp(log(1-(1/Tr))/n)
    
    Vref<-vmedia*a*(-log(1-b))^(1/k)
  }
  
  if(equation=="EWTS_G"){
    
    a<-log(n)^((1/k)-1)
    
    b<-k*gamma(1+(1/k))
    
    c<-k*log(n)-log(-log(1-(1/Tr)))
    
    Vref<-(vmedia*a*c)/b
  }
  
  if(equation=="EWTS_D"){
    
    c1<-1+((k-1)/(k*log(n)))
    
    num1<-log(k*gamma(1+(1/k))*(log(n)^((k-1)/k)))
    
    den1<-k*log(n)-(k-1)
    
    c2<-1+(num1/den1)
    
    a<-log(n)^((1/k)-1)
    
    b<-c1*k*gamma(1+(1/k))
    
    c<-c1*c2*k*log(n)
    
    d<-log(-log(1-(1/Tr)))
    
    Vref<-vmedia*(a/b)*(c-d)
  }
  
  if(equation=="Tormentas"){
    
    maximos<-Maximos_MIS(velocidad, dias, dias_sep,diezmin)
    
    maximos<-maximos[maximos>limite]
    
    ordenar<-order(maximos); maximos<-maximos[ordenar]
    
    if(length(maximos)<2){
      stop("Error: 'limite' debe ser inferior")
    }
    
    longitud<-seq(1,length(maximos),1)
    
    if(prob_method==1){
      prob<-longitud/(length(maximos)+1)
    }
    else if(prob_method==2){
      prob<-(longitud-0.44)/(length(maximos)+0.12)
    }
    else{
      stop("prob_method es un argumento inv�lido. Debe ser 1 o 2")
    }
     
    y<-(-log(-log(prob)))
    
    if(diezmin==TRUE){
      a�os<-length(velocidad)/(8760*6)
    } 
    else{
      a�os<-length(velocidad)/(8760*3600)
    }
    
    tormentas_a�o<-length(maximos)/a�os
    
    regresion<-lm(maximos~y)
    
    A<-as.vector(regresion$coefficients[1]); B<-as.vector(regresion$coefficients[2])
    
    Pref<-(Tr*tormentas_a�o)/(Tr*tormentas_a�o+1)
    
    Vref<-A+B*(-log(-log(Pref)))
    
    
    df1<-as.data.frame(matrix(nrow = length(maximos), ncol = 2)); df1$x<-y; df1$y<-maximos
    
    p6<-ggplot(data = df1, aes(x=x,y=y))+
      geom_point(size=3, colour="red", fill=I("red"))+
      geom_smooth(method = lm, formula = y~x, lwd=1.5)+
      scale_x_continuous(name = "-ln(-ln(F(x)))") +
      scale_y_continuous(name = "Wind Speed (m/s)") +
      labs(colour="Ajuste")+
      ggtitle(paste("Gumbel Fit for Independent Storms", "Vref=",round(Vref, digits = 2),"R^2=", round(summary(regresion)$r.squared, digits = 3)))  +
      theme_economist() +
      theme(legend.position = "bottom", legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.key.size = unit(1, "cm"),
            plot.title = element_text(family="Bauhaus 93", size = 18),
            text = element_text(family = "Bauhaus 93"),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 9),
            legend.title=element_text(face = "bold", size = 9, family = "Bauhaus 93"))
    
    resultado<-list(Vref=Vref, graph=p6)
    return(resultado)
    
  }
  
  if(equation=="Periodic_Maxima"){
    
    df<-as.data.frame(matrix(ncol = 2, nrow = length(velocidad))); 
    
    df$velocidad<-velocidad; df$fecha<-fecha
    
    df$fecha<-ymd(df$fecha); df$a�o<-year(df$fecha); df$mes<-month(df$fecha); 
    
    df$mes_a�o<-paste(df$mes,"-", df$a�o, sep="")
    
    if(tipo_per=="A�o"){
      
      maximos<-as.vector(tapply(df$velocidad, df$a�o, max, na.rm=TRUE ))
      
      #maximos<-maximos[maximos>limite]
      
      ordenar<-order(maximos); maximos<-maximos[ordenar]
      
      if(length(maximos)<2){
        stop("Error: 'limite' debe ser inferior")
      }
      
      longitud<-seq(1,length(maximos),1)
      
      if(prob_method==1){
        prob<-longitud/(length(maximos)+1)
      }
      else if(prob_method==2){
        prob<-(longitud-0.44)/(length(maximos)+0.12)
      }
      else{
        stop("prob_method es un argumento inv�lido. Debe ser 1 o 2")
      } 
      
      y<-(-log(-log(prob))) 
      
      if(diezmin==TRUE){
        a�os<-length(velocidad)/(8760*6)
      }
      else{
        a�os<-length(velocidad)/(8760*3600)
      }
      
      tormentas_a�o<-length(maximos)/a�os
      
      regresion<-lm(maximos~y)
      
      A<-as.vector(regresion$coefficients[1]); B<-as.vector(regresion$coefficients[2])
      
      Pref<-(Tr*tormentas_a�o)/(Tr*tormentas_a�o+1)
      
      Vref<-A+B*(-log(-log(Pref)))
    }
    
    if(tipo_per=="Mes"){
      
      maximos<-as.vector(tapply(df$velocidad, df$mes_a�o, max, na.rm=TRUE ))
      
      #maximos<-maximos[maximos>limite]
      
      ordenar<-order(maximos); maximos<-maximos[ordenar]
      
      if(length(maximos)<2){
        stop("Error: 'limite' debe ser inferior")
      }
      
      longitud<-seq(1,length(maximos),1)
      
      if(prob_method==1){
        prob<-longitud/(length(maximos)+1)
      }
      else if(prob_method==2){
        prob<-(longitud-0.44)/(length(maximos)+0.12)
      }
      else{
        stop("prob_method es un argumento inv�lido. Debe ser 1 o 2")
      } 
      
      y<-(-log(-log(prob))) 
      
      if(diezmin==TRUE){
        mes<-(length(velocidad)/(8760*6))*12
      }
      else{
        mes<-(length(velocidad)/(8760*3600))*12
      }
      
      tormentas_mes<-length(maximos)/mes
      
      regresion<-lm(maximos~y)
      
      A<-as.vector(regresion$coefficients[1]); B<-as.vector(regresion$coefficients[2])
      
      Pref<-(Tr*tormentas_mes)/(Tr*tormentas_mes+1)
      
      Vref<-A+B*(-log(-log(Pref)))
      
    }
    
    
  } 
  
  if(equation=="Gumbel"){
    
    library(ExtDist)
    
    maximos<-Maximos_MIS(velocidad, dias, dias_sep, diezmin)
    
    maximos<-maximos[maximos>limite]
    
    ordenar<-order(maximos); maximos<-maximos[ordenar]
    
    if(length(maximos)<2){
      stop("Error: 'limite' debe ser inferior")
    }
    
    longitud<-seq(1,length(maximos),1)
    
    if(prob_method==1){
      prob<-longitud/(length(maximos)+1)
    }
    else if(prob_method==2){
      prob<-(longitud-0.44)/(length(maximos)+0.12)
    }
    else{
      stop("prob_method es un argumento inv�lido. Debe ser 1 o 2")
    }
    
    y<-(-log(-log(prob)))
    
    if(diezmin==TRUE){
      a�os<-length(velocidad)/(8760*6)
    } 
    else{
      a�os<-length(velocidad)/(8760*3600)
    }
    
    tormentas_a�o<-length(maximos)/a�os
    
    regresion<-lm(maximos~y)
    
    A<-eGumbel(velocidad[is.na(velocidad)==FALSE], method = gumbel_method)
    
    scale<-A[[2]]; location<-A[[1]]
    
    Pref<-(Tr*tormentas_a�o)/(Tr*tormentas_a�o+1)
    
    Vref<-scale+location*(-log(-log(Pref)))
    
    return(Vref)
    
  }
  
  return(Vref)
}
