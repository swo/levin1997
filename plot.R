x = read_tsv('results.tsv')

p_cac = x %>%
  mutate(y=100*host) %>%
  mutate_at(vars(sS, sR), function(x) factor(round(x, 3))) %>%
  ggplot(aes(x=cac, y=y, color=factor(fnz))) +
  geom_hline(yintercept = 100, col='gray', linetype='dashed') +
  geom_point() +
  geom_line() +
  facet_grid(sS ~ sR, labeller=label_both) +
  xlab('consumption among consumers (tx per year)') +
  ylab('mean host resistance') +
  theme_classic() +
  scale_x_continuous(expand=c(0, 0), limits=c(0, 2.05)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 105)) +
  theme(aspect.ratio=0.75)

ggsave('plot_cac.pdf', plot=p_cac, useDingbats=FALSE)

p_fnz = x %>%
  mutate(y=100*host) %>%
  mutate_at(vars(sS, sR), function(x) factor(round(x, 3))) %>%
  ggplot(aes(x=fnz, y=y, color=factor(cac))) +
  geom_hline(yintercept = 100, col='gray', linetype='dashed') +
  geom_point() +
  geom_line() +
  facet_grid(sS ~ sR, labeller=label_both) +
  xlab('fraction consuming') +
  ylab('mean host resistance') +
  theme_classic() +
  #scale_x_continuous(expand=c(0, 0), limits=c(0, 2.05)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 105)) +
  theme(aspect.ratio=0.75)

ggsave('plot_fnz.pdf', plot=p_fnz, useDingbats=FALSE)
