x = read_tsv('results.tsv')

p = x %>%
  filter(s %in% c(0.005, 0.04)) %>%
  mutate(y=100*host) %>%
  ggplot(aes(x=T, y=y, group=factor(s))) +
  geom_hline(yintercept = 100, col='gray', linetype='dashed') +
  geom_point() +
  geom_line() +
  xlab('treatments per year') +
  ylab('mean resistance') +
  theme_classic() +
  scale_x_continuous(expand=c(0, 0), limits=c(0, 2.05)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 105)) +
  theme(aspect.ratio=0.75,
        text=element_text(family='Lato'))

show(p)

ggsave('plot.pdf', plot=p, useDingbats=FALSE)