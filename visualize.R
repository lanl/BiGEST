### Last updated: 01/02/2025
### Written by Taehyung Kwon & Lisa D. Adriani

# how to run this script
	# Rscript bigest_visualize.R FILE_IN PROJECT_NAME

# how to set up environment
	# conda create -n ENV_NAME conda-forge::r-base conda-forge::r-ggplot2 conda-forge::r-rcolorbrewer conda-forge::r-gggenes

initialize=function() {

	read_params=function(){

		args=commandArgs(trailingOnly=TRUE)
		
		# assign global variables
		global_vars=c(
			'file_in',
			'project'
			)
		
		for (index in 1:length(global_vars)) {
			
			# assign each global variables
			assign(global_vars[index],
				args[index],
				envir=.GlobalEnv)
			
		}

		if (is.na(file_in)) {

			print('ERROR: BGC input table not found...')
			stop()

		}
		
		if (is.na(project)) {

			project=file_in
			
		}
		
		save_dir<<-paste0(project,'/','figures/')
		
		dir.create(file.path(paste0('figures')), showWarnings=FALSE)
		
		dir.create(file.path(save_dir), showWarnings=FALSE)
		
		
		# the limit of the number of total domains
		limit_num_colors<<-20

	}

	set_params=function() {

		# lane interval between domains
		stack_interval<<-0.05

		# height_multiplier for figure height
		height_multiplier<<-3	
		
	}

	load_libraries=function() {

		suppressMessages(suppressWarnings(library(ggplot2, warn.conflicts=FALSE, quietly=TRUE)))
		suppressMessages(suppressWarnings(library(gggenes, warn.conflicts=FALSE, quietly=TRUE)))
		suppressMessages(suppressWarnings(library(RColorBrewer, warn.conflicts=FALSE, quietly=TRUE)))
		
		# set color palette
		# remove grey from set 1
		getPalette.set1<<-colorRampPalette(brewer.pal(8, 'Set1'))
		getPalette.set2<<-colorRampPalette(brewer.pal(8, 'Set2'))
		getPalette.accent<<-colorRampPalette(brewer.pal(8, 'Accent'))


	}

	# parameters
	read_params()
	set_params()


	# load libraries
	load_libraries()

}

module_assign_value=function(order.vec, value.vec) {
	
	out.vec=list()
	
	if (length(order.vec)!=length(value.vec)) {
		print('ERROR: different lengths matched')
		stop()
	}
	
	for (i in 1:length(order.vec)) {
		out.vec[[order.vec[i]]]=value.vec[i]
	}
	
	return(out.vec)

}


visualize=function() {
	
	assign_domain_class_color=function(data_in){

		# find all the domain classes
		vec_domain=na.omit(unique(as.character(data_in$domain_class)))

		# reorder vector domain
		vec_domain=vec_domain[order(vec_domain)]

		# load previously saved color list
		file_color=paste0('figures/list_domain_color.Rdata')

		# if there is previously defined color scheme, use it
		if (file.exists(file_color)){
			
			load(file_color)
			
			# if all domains are already in the list_domain_color_prev
			if (all(vec_domain %in% names(list_domain_color))) {

				
				save(list_domain_color,file=file_color)
				
				return(list_domain_color)

			} else {

				# find new domains
				vec_domain_new=vec_domain[! vec_domain %in% names(list_domain_color)]

				list_domain_color_prev=list_domain_color

				list_domain_color_prev[['-']]=NULL

			}
	
			
		} else {

			list_domain_color_prev=list()
			vec_domain_new=vec_domain
			
		}

		# exclude empty domain class ""
		vec_domain_new=vec_domain_new[which(vec_domain_new!='-')]

		# get the color vector
		vec_color=getPalette.set1(limit_num_colors)

		# pick the colors that haven't been assigned
		vec_color_not_used=vec_color[!vec_color %in% unlist(list_domain_color_prev)]

		# generate randomly sampled colors
		vec_color_new=sample(
			vec_color_not_used,
			length(vec_domain_new),
			replace=FALSE)
		
		list_domain_color_new=module_assign_value(
			c('-',vec_domain_new),
			c('grey60',vec_color_new)
			)
		
		list_domain_color=c(list_domain_color_prev,list_domain_color_new)

		save(list_domain_color,file=file_color)
		
		return(list_domain_color)

	}

	stack_domains=function(subset_in){
	
		# declare variables
		subset_in$stack=as.numeric(0)
		prev_end=0
		prev_tool=0
		max_end=0

		# reorder the dataframe
		subset_in_ordered=subset_in[order(subset_in$annotation_tool,subset_in$domain_start,subset_in$domain_length,subset_in$domain_class,subset_in$domain_name),]

		for (i in 1:nrow(subset_in_ordered)) {
			
			row=subset_in_ordered[i,]

			# for aSModule, just skip
			if ( row$module_class!='-' ) {
				
				#subset_in_ordered$stack[i]=-1*stack_interval
				subset_in_ordered$stack[i]=0
				
				next

			}

			# if this is the first domain, or the module class is not 
			if (i==1) {

				stack_lane=0
				
			} else {

				# if tool is different
				if (row$annotation_tool != prev_tool) {

					stack_lane=0
					prev_end=0
					max_end=0

				} else {

					if (row$domain_start >= max_end) {

						stack_lane=0
							
					} else if ( prev_end <= row$domain_start & row$domain_start < max_end) {

						stack_lane=stack_lane
					
					} else {
						
						stack_lane=stack_lane+stack_interval

					}

					# # if domain overlaps
					# if (row$domain_start < prev_end) {
						
					# 	stack_lane=stack_lane+stack_interval
						
					# } else if (row$domain_start >= max_end) {

					# 	stack_lane=0
							
					# }

				}
				
			}

			subset_in_ordered$stack[i]=stack_lane
			max_end=max(row$domain_end,prev_end)
			prev_end=row$domain_end
			prev_tool=row$annotation_tool
			
		}
	
		return(subset_in_ordered)

	}

	visualize_each_cluster=function(bgc_id,subset_in,list_domain_color) {
		
		file_fig=paste0(save_dir,'/',bgc_id,'.pdf')

		subset_stack=stack_domains(subset_in)
		
		figure_height=0

		# estimate figure size
		for (annotation_tool in levels(subset_stack$annotation_tool)) {
			
			# total number of stacks
			num_stack=length(unique(subset_stack$stack[which(subset_stack$annotation_tool==annotation_tool)]))

			# add 2 to the total number (additional figure features)
			figure_height=figure_height+(num_stack+2)*height_multiplier

		}

		fig_main=ggplot(subset_stack)+
			
			geom_gene_arrow(
				aes(
					xmin=domain_start,
					xmax=domain_end,
					y=stack,
					fill=domain_class,
					color=module_class,
					forward=orientation
					),
				
				arrow_body_height=unit(2.5, "mm"),
				arrowhead_height=unit(2.5, "mm"),
				arrowhead_width=unit(1, "mm"),
				alpha=0.8
				) +

			geom_text(
				aes(
					x=domain_start,
					y=stack+stack_interval*0.5,
					label=domain_name,
					),
					vjust=0.5,
					hjust=0,
					size=1
				
				)+
			
			scale_fill_manual(
				values=list_domain_color,
				name='',
				na.value=NA,
				drop=TRUE
				)+

			scale_color_manual(
				values=c(
					'PKS Module'='firebrick1',
					'NRPS Module'='dodgerblue',
					'-'=NA
					),
					name='',na.value=NA,drop=TRUE
				)+
			
			scale_y_continuous(expand=c(0,stack_interval*0.6))+
			scale_x_continuous(breaks=seq(
				floor(min(c(subset_stack$domain_start,subset_stack$domain_end))/10000)*10000,
				ceiling(max(c(subset_stack$domain_start,subset_stack$domain_end))/10000)*10000,
				10000
				),
				expand=c(0.02,0)
				)+
			
			labs(x='genomic location',y='',title=bgc_id)+
			
			guides(
				fill='none',
				color='none',
				alpha='none')+

			facet_grid(annotation_tool~.,scales='free_y',space='free_y')+
			
			theme(
				text=element_text(size=5,hjust=0.5,vjust=0.5),
				line=element_line(linewidth=0.5),
				panel.background=element_blank(),
				panel.border=element_rect(linewidth=0.5,color='black',linetype='solid',fill=NA),
				panel.grid.major.x=element_line(linewidth=0.25,color='grey66',linetype='dashed'),
				panel.spacing=unit(1,'mm'),
				plot.title=element_text(size=5),
				axis.ticks.y=element_blank(),
				axis.line.x=element_line(linewidth=0.5),
				axis.text.y=element_blank(),
				axis.title=element_blank(),
				plot.margin=unit(c(1,1,1,1),'mm'),
				strip.background=element_rect(color=NA,linewidth=0.5,fill='white'),
				strip.text=element_text(
					colour="black",size=3,face='bold',margin=margin(1,1,1,1,'mm'),angle=0,hjust=.5,vjust=.5),
				legend.background=element_rect(color=NA,linewidth=0.3),
				legend.key=element_rect(fill='NA'),
				legend.key.size=unit(0.2,'lines'),
				legend.position='bottom',
				legend.box='vertical')

		ggsave(file_fig, fig_main, unit='mm', width=170, height=figure_height,dpi=600)

	}
	
	modify_domain_name=function(domain_name_vec) {
		
		# remove PKS or NRPS from domain names
		domain_name_vec=gsub('PKS_','',as.character(domain_name_vec))

		domain_name_vec=gsub('_NRPS','',as.character(domain_name_vec))

		return(domain_name_vec)

	}

	modify_domain_class=function(data_in) {

		# add synthaser domain class to matching antiSMASH domains
		df_class_to_name=data_in[which(data_in$annotation_tool=='BiGEST'),c('domain_name','domain_class')]
		colnames(df_class_to_name)=c('domain_name','domain_class_BiGEST')

		df_class_to_name_dedup=df_class_to_name[!duplicated(df_class_to_name),]
		
		data_merge=merge(data_in,df_class_to_name_dedup,by='domain_name',all.x = TRUE)

		# replace domain class for domains (except for modules)
		data_merge$domain_class[which(data_merge$domain_name!='-')]=data_merge$domain_class_BiGEST[which(data_merge$domain_name!='-')]

		# factorize domain class
		data_merge$domain_class=as.factor(data_merge$domain_class)

		return(data_merge)

	}

	load_data=function() {

		# load BiGEST output
		data_raw=data.frame(
			read.table(
				file_in,
				header=F,
				sep='\t',
				quote=NULL
				)
			)

		# assign column names
		colnames(data_raw)=c('bgc_id','annotation_tool','domain_class','domain_name','module_class','domain_start','domain_end','domain_strand')

		### temp: remove deepBGC features
		data_raw=data_raw[which(data_raw$annotation_tool != 'deepbgc'),]

		data_raw$domain_class=as.character(data_raw$domain_class)
		
		# change empty domain_class values to '-'
		data_raw$domain_class[which(data_raw$domain_class=='' | is.na(data_raw$domain_class))]='-'
		
		data_raw$domain_class[which(data_raw$module_class!='')]=NA

		# change empty module_class values to '-'
		data_raw$module_class[which(data_raw$module_class=='' | is.na(data_raw$module_class))]='-'
		
		# add orientation column
		data_raw$orientation=ifelse(data_raw$domain_strand=='+',1,0)

		# factorize domain type
		data_raw$annotation_tool=factor(data_raw$annotation_tool,levels=c('BiGEST','antiSMASH','gecco','deepbgc'))

		# domain length (zero-based)
		data_raw$domain_length=abs(data_raw$domain_end-data_raw$domain_start)

		# modify domain names
		data_raw$domain_name=modify_domain_name(data_raw$domain_name)
		
		### temp
		# # modify domain classes
		# data_class_modified=modify_domain_class(data_raw)

		return(data_raw)
	}

	# main function
	# call data file
	data_in=load_data()

	# set domain color
	list_domain_color=assign_domain_class_color(data_in)

	for (bgc_each in unique(data_in$bgc_id)) {

		print(paste0('INFO: visualize | BGC ID=',bgc_each))
		# subset dataframe for each BGC
		subset_in=data_in[which(data_in$bgc_id==bgc_each),]

		# visualize each BGC
		visualize_each_cluster(bgc_each,subset_in,list_domain_color)

	}

}

initialize()
visualize()