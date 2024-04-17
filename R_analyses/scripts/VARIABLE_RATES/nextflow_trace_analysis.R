library(tidyverse)

trace <- read.table("~/Documents/PhD/Projects/2023/Avian_scRNAseq/Meetings/Alison07112023/trace.txt", fill = T)
trace <- read.table(args[1])

colnames(trace) <- trace[1,]
colnames(trace)[colnames(trace) == ""] <- c(1:100)
trace <- trace[-1,]

remove_Ns <- trace[trace$name == "remove_Ns",]
remove_Ns %>% 
  group_by(exit) %>% 
  summarise(n = n())

trace_sum <- trace %>% 
  group_by(name, exit) %>% 
  summarise(n = n())
trace_sum <- trace_sum[trace_sum$exit != "",]

order_process <- unique(trace$name)
order_process <- order_process[order_process != ""]

#order trace_sum df by the name but using order of order_processes
trace_sum <- trace_sum[c(9,3,4,8,7,2,10,11,15,16,12,6,17,13,14,5,1),]
trace_sum$exit[trace_sum$exit %in% c(0, "CACHED", "COMPLETED")] <- "PASS" 
trace_sum$nname <- paste(letters[1:17], trace_sum$name)

#order ggplot bar in order of appearance in df
trace_sum_plot <- trace_sum %>% ggplot(aes(x = nname, y = n, fill = exit)) + 
  geom_bar(stat = 'identity')  +# add number on top of each bar
  geom_text(aes(label = n), vjust = -0.5, size = 3) 

#save plot with ggsave 
ggsave("plots/trace_sum_plot.pdf", plot = trace_sum_plot, width = 30, height = 5, units = "in", dpi = 300)




