library(decentralized1bit)

library(ggplot2)
library(ggthemes)

error_data_positive <- data.frame(
  Local = positive[,2],  
  Avg = positive[,3],       
  subGD = positive[,4],  
  Our = positive[,5]    
)

plot_data1 <- reshape2::melt(error_data_positive, 
                            variable.name = "Method",
                            value.name = "Error")

ggplot(plot_data1, aes(x = Method, y = Error, fill = Method)) + 
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() + 
  theme(legend.position = "none", axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24)) + 
  labs(x = "", y = "") +
  scale_y_continuous(limits = c(0,1)) 
 
  

##################################
error_data_neutral <- data.frame(
  Local = neutral[,2],  
  Avg = neutral[,3],       
  subGD = neutral[,4],  
  Our = neutral[,5]    
)

plot_data2 <- reshape2::melt(error_data_neutral, 
                             variable.name = "Method",
                             value.name = "Error")

ggplot(plot_data2, aes(x = Method, y = Error, fill = Method)) + 
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() + 
  theme(legend.position = "none", axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))  + 
  labs(x = "", y = "") +
  scale_y_continuous(limits = c(0,1))



##################################
error_data_negative <- data.frame(
  Local = negative[,2],  
  Avg = negative[,3],       
  subGD = negative[,4],  
  Our = negative[,5]    
)

plot_data3 <- reshape2::melt(error_data_negative, 
                             variable.name = "Method",
                             value.name = "Error")

ggplot(plot_data3, aes(x = Method, y = Error, fill = Method)) + 
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() + 
  theme(legend.position = "none", axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))  + 
  labs(x = "", y = "") +
  scale_y_continuous(limits = c(0,1))




#############################################
error_data_eeg1 <- data.frame(
  Local = eeg1[,2],  
  Avg = eeg1[,3],       
  subGD = eeg1[,4],  
  Our = eeg1[,5]    
)

plot_data4 <- reshape2::melt(error_data_eeg1, 
                             variable.name = "Method",
                             value.name = "Error")

ggplot(plot_data4, aes(x = Method, y = Error, fill = Method)) + 
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() + 
  theme(legend.position = "none", axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))  + 
  labs(x = "", y = "") +
  scale_y_continuous(limits = c(0,1))



#############################################
error_data_eeg2 <- data.frame(
  Local = eeg2[,2],  
  Avg = eeg2[,3],       
  subGD = eeg2[,4],  
  Our = eeg2[,5]    
)

plot_data5 <- reshape2::melt(error_data_eeg2, 
                             variable.name = "Method",
                             value.name = "Error")

ggplot(plot_data5, aes(x = Method, y = Error, fill = Method)) + 
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() + 
  theme(legend.position = "none", axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))  + 
  labs(x = "", y = "") +
  scale_y_continuous(limits = c(0,1))


#############################################
error_data_eeg3 <- data.frame(
  Local = eeg3[,2],  
  Avg = eeg3[,3],       
  subGD = eeg3[,4],  
  Our = eeg3[,5]    
)

plot_data6 <- reshape2::melt(error_data_eeg3, 
                             variable.name = "Method",
                             value.name = "Error")

ggplot(plot_data6, aes(x = Method, y = Error, fill = Method)) + 
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() + 
  theme(legend.position = "none", axis.text.x = element_text(size = 24), axis.text.y = element_text(size = 24))  + 
  labs(x = "", y = "") +
  scale_y_continuous(limits = c(0,1))

