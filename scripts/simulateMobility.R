library(plyr)

code = '
Rcpp::List simulate_fun_separate(int nday, arma::vec onsetstates, arma::vec onsetweeks, double r0, arma::mat mobility_commuting, arma::mat mobility_airtravel, arma::vec pops, double mob_rate) {
  
  onsetweeks = onsetweeks - min(onsetweeks);
  int ncluster = onsetweeks.n_elem;
  
  arma::mat I_arr_out_commute = arma::mat(49*ncluster*nday,49);
  arma::mat I_arr_out_air = arma::mat(49*ncluster*nday,49);

  arma::mat inc_arr = arma::mat(49*ncluster, nday).fill(0.);
  
  arma::mat I_arr_commute = arma::mat(49*ncluster, 49).fill(0.);
  arma::mat I_arr_air = arma::mat(49*ncluster, 49).fill(0.);
  arma::mat I_mat_local = arma::mat(ncluster,49).fill(0.);
  
  arma::mat S_mat_commute = arma::mat(49,49).fill(0.);
  arma::mat S_mat_air = arma::mat(49,49).fill(0.);
  arma::rowvec S_vec_local = arma::rowvec(49).fill(0.);
  for (int i = 0; i < 49; i++){
    S_vec_local(i) = pops(i);
  }
  
  arma::mat N_mat_commute = arma::mat(49,49).fill(0.);
  arma::mat N_mat_air = arma::mat(49,49).fill(0.);
  arma::rowvec N_vec_local = arma::rowvec(49).fill(0.);
  for (int i = 0; i < 49; i++){
    N_vec_local(i) = pops(i);
  }

  double gamma = 0.25;
  double beta = gamma*r0;

  double return_rate = 1;
  
  arma::mat mob_mat_commute = mobility_commuting * mob_rate;
  arma::mat mob_mat_air = mobility_airtravel * mob_rate;

  arma::mat tot_infections_commute =  arma::mat(49, 49);
  arma::mat tot_infections_air =  arma::mat(49, 49);
  arma::rowvec tot_infections_local =  arma::rowvec(49);
  
  arma::rowvec I_vec_local =  arma::rowvec(49);
  arma::mat I_mat_commute =  arma::mat(49, 49);
  arma::mat I_mat_air =  arma::mat(49, 49);
  arma::mat I_leave_commute =  arma::mat(49, 49);
  arma::mat I_leave_air =  arma::mat(49, 49);
  arma::mat I_return_commute =  arma::mat(49, 49);
  arma::mat I_return_air =  arma::mat(49, 49);
  arma::mat mobility_I_commute =  arma::mat(49, 49);
  arma::mat mobility_I_air =  arma::mat(49, 49);

  arma::mat mobility_I =  arma::mat(49, 49);
  arma::vec div = arma::vec(49);
  
  arma::mat infections_commute = arma::mat(49,49);
  arma::mat infections_air = arma::mat(49,49);
  arma::rowvec infections_local = arma::rowvec(49);

  arma::mat recoveries_commute = arma::mat(49,49);
  arma::mat recoveries_air = arma::mat(49,49);
  arma::rowvec recoveries_local = arma::rowvec(49);

  arma::mat S_leave_commute =  arma::mat(49, 49);
  arma::mat S_return_commute =  arma::mat(49, 49);
  arma::mat S_leave_air =  arma::mat(49, 49);
  arma::mat S_return_air =  arma::mat(49, 49);
  arma::mat mobility_S_commute =  arma::mat(49, 49);
  arma::mat mobility_S_air =  arma::mat(49, 49);
  
  arma::mat N_leave_commute =  arma::mat(49, 49);
  arma::mat N_return_commute =  arma::mat(49, 49);
  arma::mat N_leave_air =  arma::mat(49, 49);
  arma::mat N_return_air =  arma::mat(49, 49);
  arma::mat mobility_N =  arma::mat(49, 49);
  arma::mat mobility_N_commute =  arma::mat(49, 49);
  arma::mat mobility_N_air =  arma::mat(49, 49);
  
  arma::rowvec save_col_commute = arma::rowvec(49);
  arma::rowvec save_col_air= arma::rowvec(49);
  arma::rowvec save_row_commute = arma::rowvec(49);
  arma::rowvec save_row_air = arma::rowvec(49);

  for (int i = 0; i < nday; i++){
    for (int j = 0; j < ncluster; j++){
    
      if (onsetweeks(j) == i){
      
        I_mat_local(j,onsetstates(j)) = pops(onsetstates(j))*0.00001;
      
      }
    }

    tot_infections_local.zeros();
    
    for (int cl = 0; cl < ncluster; cl++){

      I_mat_commute = I_arr_commute.submat(cl*49,0,(cl+1)*49-1,48);
      I_mat_air = I_arr_air.submat(cl*49,0,(cl+1)*49-1,48);
      I_vec_local = I_mat_local.row(cl);
      div = (arma::sum(I_mat_commute,1) + I_vec_local.t() + arma::sum(I_mat_air,1)) / (arma::sum(N_mat_commute,1) + arma::sum(N_mat_air,1) + N_vec_local.t());

      for (int j = 0; j < 49; j++){
        
        infections_commute.col(j) = beta * S_mat_commute.col(j) * div(j);
        infections_air.col(j) = beta * S_mat_air.col(j) * div(j);
        infections_local(j) = beta * S_vec_local(j) * div(j);

        I_leave_commute.col(j) = mob_mat_commute.col(j) * I_vec_local(j);
        I_leave_air.col(j) = mob_mat_air.col(j) * I_vec_local(j);

      }

      tot_infections_commute += infections_commute;
      tot_infections_air += infections_air;
      tot_infections_local += infections_local;
      
      I_return_commute = return_rate * I_mat_commute;
      I_return_air = return_rate * I_mat_air;
      
      mobility_I_commute = I_leave_commute - I_return_commute;
      mobility_I_air = I_leave_air - I_return_air;

      save_col_commute = -1.0 * arma::sum(I_leave_commute,0) + arma::sum(I_return_commute,0);
      save_col_air = -1.0 * arma::sum(I_leave_air,0) + arma::sum(I_return_air,0);
      
      recoveries_commute = gamma * I_mat_commute;
      recoveries_air = gamma * I_mat_air;
      recoveries_local = gamma * I_vec_local;

      I_mat_commute += infections_commute + mobility_I_commute - recoveries_commute;
      I_mat_air += infections_air + mobility_I_air - recoveries_air;
      I_vec_local += infections_local + save_col_commute + save_col_air - recoveries_local;

      I_mat_commute.elem(find(I_mat_commute<0)).fill(0.0);
      I_mat_air.elem(find(I_mat_air<0)).fill(0.0);
      I_vec_local.elem(find(I_vec_local<0)).fill(0.0);

      I_arr_commute.submat(cl*49,0,(cl+1)*49-1,48) = I_mat_commute;
      I_arr_air.submat(cl*49,0,(cl+1)*49-1,48) = I_mat_air;
      I_mat_local.row(cl) = I_vec_local;

      inc_arr.submat(cl*49,i,(cl+1)*49-1,i) = infections_local.t();
      
      I_arr_out_commute.submat((cl*49+i*49*ncluster),0,(cl+1)*49-1+i*49*ncluster,48) = I_mat_commute;
      I_arr_out_air.submat((cl*49+i*49*ncluster),0,(cl+1)*49-1+i*49*ncluster,48) = I_mat_air;

    }
    
    for (int j = 0; j < 49; j++){
      
      S_leave_commute.col(j) = mob_mat_commute.col(j) * (double) S_vec_local(j);
      S_leave_air.col(j) = mob_mat_air.col(j) * (double) S_vec_local(j);

      N_leave_commute.col(j) = mob_mat_commute.col(j) * N_vec_local(j);
      N_leave_air.col(j) = mob_mat_air.col(j) * N_vec_local(j);

    }
    

    S_return_commute = return_rate * S_mat_commute;
    S_return_air = return_rate * S_mat_air;

    mobility_S_commute = S_leave_commute - S_return_commute;
    mobility_S_air = S_leave_air - S_return_air;

    save_col_commute = -1.0 * arma::sum(S_leave_commute,0) + arma::sum(S_return_commute,0);
    save_col_air = -1.0 * arma::sum(S_leave_air,0) + arma::sum(S_return_air,0);

    for (int i = 0; i < 49; i++){
      S_vec_local(i) = S_vec_local(i) - tot_infections_local(i) + save_col_commute(i) + save_col_air(i);
    }
    
    S_mat_commute = S_mat_commute - tot_infections_commute + mobility_S_commute;
    S_mat_air = S_mat_air - tot_infections_air + mobility_S_air;

      S_mat_commute.elem(find(S_mat_commute<0)).fill(0.0);
      S_mat_air.elem(find(S_mat_air<0)).fill(0.0);
      S_vec_local.elem(find(S_vec_local<0)).fill(0.0);

    
    N_return_commute = return_rate * N_mat_commute;
    N_return_air = return_rate * N_mat_air;

    mobility_N_commute = N_leave_commute - N_return_commute;
    mobility_N_air = N_leave_air - N_return_air;

    save_col_commute = -1.0 * arma::sum(N_leave_commute,0) + arma::sum(N_return_commute,0);
    save_col_air = -1.0 * arma::sum(N_leave_air,0) + arma::sum(N_return_air,0);

    for (int i = 0; i < 49; i++){
      N_vec_local(i) = N_vec_local(i) + save_col_commute(i) + save_col_air(i);
    }
    
    N_mat_commute = N_mat_commute - tot_infections_commute + mobility_N_commute;
    N_mat_air = N_mat_air - tot_infections_air + mobility_N_air;

  }
  return(Rcpp::List::create(Rcpp::Named("inc")=inc_arr,Rcpp::Named("commuting")=I_arr_out_commute,Rcpp::Named("air")=I_arr_out_air));
}

'

Rcpp::cppFunction(code=code, env=.GlobalEnv, depends="RcppArmadillo")


getMobilityMatrix <- function(){
  
  states = c(state.name[-c(2,11)],"District Of Columbia")
  states_abb = c(state.abb[match(state.name[-c(2,11)],state.name)],"DC")
  
  commuting_data = read.csv("~/Downloads/commuting_flows.csv",h=F)
  colnames(commuting_data) = c("residence.state.code","residence.county.code","residence.state","residence.county",
                               "workplace.state.code","workplace.county.code","workplace.state","workplace.county",
                               "commuting.flow","commuting.flow.moe")
  
  commuting_data = commuting_data[,c(3,7,9)]
  colnames(commuting_data) = c("ORIG","DEST","COUNT")
  commuting_data[commuting_data == "District of Columbia"] = "District Of Columbia"
  commuting_data = commuting_data[commuting_data$ORIG %in% states & commuting_data$DEST %in% states,]
  commuting_data = plyr::ddply(commuting_data,.(ORIG,DEST),summarise,N = sum(COUNT))
  commuting_data = commuting_data[2:nrow(commuting_data),]
  commuting_data = acast(commuting_data,ORIG~DEST,value.var="N")
  commuting_data[is.na(commuting_data)] = 1
  commuting_data[commuting_data==0] = 1
  
  commuting_data = (commuting_data+t(commuting_data))/2
  
  for (i in ncol(commuting_data)){
    commuting_data[i,i] = NA
  }
  commuting_data = commuting_data[match(states,colnames(commuting_data)),match(states,colnames(commuting_data))]
  
  airtravel_data = read.csv("~/Downloads/T_T100D_MARKET_ALL_CARRIER-2.csv")
  airtravel_data = airtravel_data[airtravel_data$ORIGIN_STATE_ABR %in% states_abb & airtravel_data$DEST_STATE_ABR %in% states_abb,]
  airtravel_data = ddply(airtravel_data,.(ORIGIN_STATE_ABR,DEST_STATE_ABR),summarise,N = sum(PASSENGERS))
  airtravel_data = acast(airtravel_data,ORIGIN_STATE_ABR~DEST_STATE_ABR,value.var="N")
  airtravel_data[is.na(airtravel_data)] = 1
  airtravel_data[airtravel_data==0] = 1
  for (i in 1:ncol(airtravel_data)){
    airtravel_data[i,i] = 0
  }
  airtravel_data = (airtravel_data+t(airtravel_data))/2
  airtravel_data = rbind(airtravel_data,rep(1,48))
  airtravel_data = cbind(airtravel_data,rep(1,49))
  colnames(airtravel_data)[49] = "DC"
  
  idxes = match(states_abb,colnames(airtravel_data))
  airtravel_data = airtravel_data[idxes,]
  airtravel_data = airtravel_data[,idxes]
  for (i in 1:ncol(airtravel_data)){
    airtravel_data[i,i] = NA
  }
  
  d = commuting_data 
  for (i in 1:ncol(d)){d[i,i] = 0}
  d2_commuting = d
  
  d = (airtravel_data/365) 
  for (i in 1:ncol(d)){d[i,i] = 0}
  d2_airtravel= d
  
  d = pmax(d2_commuting,d2_airtravel)
  for (i in 1:ncol(d)){d[,i] = d[,i]/pops_raw[i]}
  
  for (i in 1:ncol(d2_commuting)){d2_commuting[,i] = d2_commuting[,i]/pops_raw[i]}
  for (i in 1:ncol(d2_airtravel)){d2_airtravel[,i] = d2_airtravel[,i]/pops_raw[i]}
  
  adjacency_mat = read.csv("/Users/simondejong/usa/adjacency_matrix.csv")
  adjacency_mat = adjacency_mat[,-c(2,11)]
  adjacency_mat = adjacency_mat[-c(2,11),]
  colnames(adjacency_mat) = states[1:48]
  rownames(adjacency_mat) = states[1:48]
  adjacency_mat = rbind(adjacency_mat,0)
  adjacency_mat = cbind(adjacency_mat,0)
  colnames(adjacency_mat)[49] = "District Of Columbia"
  rownames(adjacency_mat)[49] = "District Of Columbia"
  adjacency_mat[49,match(c("Virginia","Maryland"),colnames(adjacency_mat))] = 1
  adjacency_mat[match(c("Virginia","Maryland"),colnames(adjacency_mat)),49] = 1
  
  d_commuting = d
  d_commuting[adjacency_mat==F] = 0
  d_airtravel = d
  d_airtravel[adjacency_mat==T] = 0
  
  return(list(d_airtravel,d_commuting,d2_airtravel,d2_commuting,adjacency_mat))
}


getPcts <- function(sim){
  states = c(state.name[-c(2,11)],"District Of Columbia")
  ncluster = dim(sim[[1]])[1]/49
  nstate = 49
  nday = dim(sim[[1]])[2]
  simout1 = array(0,dim=c(nstate,nday,ncluster))
  for (clidx in 1:ncluster){
    mat = sim$inc[((clidx-1)*nstate+1):(clidx*nstate),]
    simout1[,,clidx] = mat
  }
  sim1 = simout1
  all_pct = matrix(0,length(states),ncluster)
  all_ons = matrix(0,length(states),ncluster)
  rs = rowSums(simout1)
  for (i in 1:ncluster){
    out_inc = simout1[,,i]
    ons = unlist(lapply(1:nstate,function(x)which(cumsum(out_inc[x,])>((rs[x])*0.05))[1]))/7
    pct = rowSums(out_inc)/rowSums(sim1)
    all_pct[,i] = pct
    all_ons[,i] = ons
  }
  return(list(all_pct,all_ons))
}



nday = 365
r0 = 1.35

states = c(state.name,"District Of Columbia")[-c(2,11)]
statepops = read.csv("/Users/simondejong/usa/statepops.csv")
statepops = statepops[statepops$year=="2012",]
statepops = statepops[statepops$ages=="total",]
pops = statepops[statepops$state.region %in% c(state.abb),]$population
pops = c(pops,statepops[statepops$state.region == "DC",]$population)
pops = pops[-c(2,11)]
pops_raw = pops
pops = pops/mean(pops)

mob = getMobilityMatrix()
d_airtravel = mob[[1]]
d_commuting = mob[[2]]



