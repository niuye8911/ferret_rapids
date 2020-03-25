/* AUTORIGHTS
Copyright (C) 2007 Princeton University

This file is part of Ferret Toolkit.

Ferret Toolkit is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#include "rsdgMissionWrapper.h"
#include <../image/image.h>
#include <cass.h>
#include <cass_stat.h>
#include <cass_timer.h>
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

#define MAXR 100
#define IMAGE_DIM 14

const char *db_dir = NULL;
const char *table_name = NULL;
const char *query_dir = NULL;
const char *output_path = NULL;
int RSDG = 0;
int totSec = 0;
FILE *fout;
int CONT = 1;
char *config_file;
// create a new rsdgMission
void *ferretMission;
void *iterationPara;
void *hashPara;
void *probePara;
void *hashContPara;
void *probeContPara;
void *iterationContPara;

int top_K = 10;

char *extra_params = "-L 8 - T 20";

int curHashNum = 8;
int curProbeNum = 20;
int curItr = 500;

// char *extra_params = "-L 10 -T 40";

int input_end;
pthread_cond_t done;
pthread_mutex_t done_mutex;

cass_env_t *env;
cass_table_t *table;
cass_table_t *query_table;

int vec_dist_id = 0;
int vecset_dist_id = 0;

int NTHREAD = 1;
int DEPTH = 1;

/* ------- The Helper Functions ------- */
char path[BUFSIZ];

// RSDG:modified signature to add configurable variable
void do_query(const char *);
int scan_dir(const char *, char *head);

// RSDG: a helper timer to measure the cost
long long getCurrentMilli() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  long long mslong = (long long)(tp.tv_sec * 1000 + tp.tv_usec / 1000);
  return mslong;
}

int dir_helper(char *dir, char *head) {
  DIR *pd = NULL;
  struct dirent *ent = NULL;
  int result = 0;
  pd = opendir(dir);
  if (pd == NULL)
    goto except;
  // first count the total number of files
  int totalCount = 0;
  while ((ent = readdir(pd)) != NULL) {
    if (ent->d_name[0] != '.') { /* If the entry is a regular file */
      totalCount++;
    }
  }
  printf("Total number of Files:%d\n", totalCount);
  if (RSDG == 1)
    setUnit(ferretMission, totalCount);
  closedir(pd);
  pd = opendir(dir);
  // setup a counter of how many iterations have been finished
  if (RSDG == 1)
    start(ferretMission);

  for (;;) {
    ent = readdir(pd);
    if (ent == NULL) { // the end
      break;
    }
    if (ent && ent->d_name[0] == '.') {
      continue;
    }
    if (scan_dir(ent->d_name, head) != 0)
      return -1;
    if (RSDG == 1) {
      finish_one_unit(ferretMission);
    }
  }
  if (RSDG == 1) {
    finish(ferretMission);
  }
  goto final;

except:
  result = -1;
  perror("Error:");
final:
  if (pd != NULL)
    closedir(pd);
  return result;
}

int scan_dir(const char *dir, char *head) {
  struct stat st;
  int ret;
  /* test for . and .. */
  if (dir[0] == '.') {
    if (dir[1] == 0)
      return 0;
    else if (dir[1] == '.') {
      if (dir[2] == 0)
        return 0;
    }
  }

  /* append the name to the path */
  strcat(head, dir);
  ret = stat(path, &st);
  if (ret != 0) {
    perror("Error:");
    return -1;
  }
  if (S_ISREG(st.st_mode))
    do_query(path);
  else if (S_ISDIR(st.st_mode)) {
    strcat(head, "/");
    dir_helper(path, head + strlen(head));
  }

  /* removed the appended part */
  head[0] = 0;
  return 0;
}

/* ------ The Stages ------ */
void scan(void) {
  const char *dir = query_dir;

  path[0] = 0;

  if (strcmp(dir, ".") == 0) {
    dir_helper(".", path);
  } else {
    scan_dir(dir, path);
  }
}

void do_query(const char *name) {
  cass_dataset_t ds;
  cass_query_t query;
  cass_result_t result;
  cass_result_t *candidate;

  // RSDG: the step below is not compuational expensive
  {

    unsigned char *HSV, *RGB;
    unsigned char *mask;
    int width, height, nrgn;
    int r;

    r = image_read_rgb_hsv(name, &width, &height, &RGB, &HSV);
    assert(r == 0);

    image_segment(&mask, &nrgn, RGB, width, height);

    image_extract_helper(HSV, mask, width, height, nrgn, &ds);

    /* free image & map */
    free(HSV);
    free(RGB);
    free(mask);
  }
  // RSDG: setting up the query
  memset(&query, 0, sizeof query);
  query.flags = CASS_RESULT_LISTS | CASS_RESULT_USERMEM;

  query.dataset = &ds;
  query.vecset_id = 0;

  query.vec_dist_id = vec_dist_id;

  query.vecset_dist_id = vecset_dist_id;

  query.topk = 2 * top_K;

  // RSDG: modified by RSDG to suppport different para
  char *extra_params_rsdg = (char *)malloc(20 * sizeof(char));
  sprintf(extra_params_rsdg, "-L %d -T %d\0", curHashNum, curProbeNum);

  query.extra_params = &extra_params_rsdg[0];

  cass_result_alloc_list(&result, ds.vecset[0].num_regions, query.topk);

  // this query get the candidate
  cass_table_query(table, &query, &result, curItr);

  memset(&query, 0, sizeof query);

  query.flags = CASS_RESULT_LIST | CASS_RESULT_USERMEM | CASS_RESULT_SORT;
  query.dataset = &ds;
  query.vecset_id = 0;

  query.vec_dist_id = vec_dist_id;

  query.vecset_dist_id = vecset_dist_id;

  query.topk = top_K;

  query.extra_params = NULL;

  // RSDG: getting the top-2K candidates(IMPORTANT)
  candidate = cass_result_merge_lists(
      &result, (cass_dataset_t *)query_table->__private, 0);
  query.candidate = candidate;

  cass_result_free(&result);

  cass_result_alloc_list(&result, 0, top_K);
  // RSDG: query the table with the candidates, with last arg being the maximum
  // iteration for computing emd
  cass_table_query(query_table, &query, &result, curItr);

  cass_result_free(candidate);
  free(candidate);
  cass_dataset_release(&ds);

  fprintf(fout, "%s", name);
  // RSDG: print result to output
  ARRAY_BEGIN_FOREACH(result.u.list, cass_list_entry_t p) {
    char *obj = NULL;
    if (p.dist == MAXDOUBLE)
      continue;
    cass_map_id_to_dataobj(query_table->map, p.id, &obj);
    assert(obj != NULL);
    fprintf(fout, "\t%s:%g", obj, p.dist);
  }
  ARRAY_END_FOREACH;

  fprintf(fout, "\n");

  cass_result_free(&result);
}

void *change_Iteration_Num(void *arg);
void *change_Hash_Num(void *arg);
void *change_Probe_Num(void *arg);
void setupMission();

int main(int argc, char *argv[]) {
  stimer_t tmr;
  int ret, i;

#ifdef PARSEC_VERSION
#define __PARSEC_STRING(x) #x
#define __PARSEC_XSTRING(x) __PARSEC_STRING(x)
  printf(
      "PARSEC Benchmark Suite Version "__PARSEC_XSTRING(PARSEC_VERSION) "\n");
  fflush(NULL);
#else
  printf("PARSEC Benchmark Suite\n");
  fflush(NULL);
#endif // PARSEC_VERSION
#ifdef ENABLE_PARSEC_HOOKS
  __parsec_bench_begin(__parsec_ferret);
#endif

  if (argc < 8) {
    printf("%s <database> <table> <query dir> <top K> <depth> <n> <out>\n",
           argv[0]);
    return 0;
  }

  db_dir = argv[1];
  table_name = argv[2];
  query_dir = argv[3];
  top_K = atoi(argv[4]);

  DEPTH = atoi(argv[5]);
  NTHREAD = atoi(argv[6]);
  if (NTHREAD != 1) {
    printf("n must be 1 (serial version)\n");
    exit(1);
  }

  output_path = argv[7];

  if (argc > 8) {
    fprintf(stdout, "RAPID is involved\n");
    // RAPID is involved
    // argv[8] = -rsdg, argv[9] = -b, argv[10] = budget
    // argv[11] = -l, argv[12] = curHashNum
    // argv[13] = -t, argv[14] = curProbeNum
    // argv[15] = -itr, argv[16] = curItr
    int cur_argc = 8;
    while (cur_argc < argc) {
      if (!strcmp(argv[cur_argc], "-rsdg")) {
        RSDG = 1;
        config_file = argv[++cur_argc];
        cur_argc++;
      } else if (!strcmp(argv[cur_argc], "-l")) {
        curHashNum = atoi(argv[++cur_argc]);
        cur_argc++;
      } else if (!strcmp(argv[cur_argc], "-t")) {
        curProbeNum = atoi(argv[++cur_argc]);
        cur_argc++;
      } else if (!strcmp(argv[cur_argc], "-itr")) {
        curItr = atoi(argv[++cur_argc]);
        cur_argc++;
      }
    }
  }

  fout = fopen(output_path, "w");
  assert(fout != NULL);

  cass_init();

  ret = cass_env_open(&env, db_dir, 0);
  if (ret != 0) {
    printf("ERROR: %s\n", cass_strerror(ret));
    return 0;
  }

  vec_dist_id = cass_reg_lookup(&env->vec_dist, "L2_float");
  assert(vec_dist_id >= 0);

  vecset_dist_id = cass_reg_lookup(&env->vecset_dist, "emd");
  assert(vecset_dist_id >= 0);

  i = cass_reg_lookup(&env->table, table_name);

  table = query_table = cass_reg_get(&env->table, i);

  i = table->parent_id;

  if (i >= 0) {
    query_table = cass_reg_get(&env->table, i);
  }

  if (query_table != table)
    cass_table_load(query_table);
  cass_map_load(query_table->map);
  cass_table_load(table);

  image_init(argv[0]);

  stimer_tick(&tmr);

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_begins();
#endif
  // setup RAPID before executing the work
  if (RSDG == 1) {
    setupMission();
  }
  scan();
#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_ends();
#endif

  stimer_tuck(&tmr, "QUERY TIME");

  ret = cass_env_close(env, 0);
  if (ret != 0) {
    printf("ERROR: %s\n", cass_strerror(ret));
    return 0;
  }

  cass_cleanup();

  image_cleanup();

  fclose(fout);

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_bench_end();
#endif

  return 0;
}

void *change_Hash_Num(void *arg) {
  int hashNum = getParaVal(hashPara);
  int newHashNum;
  newHashNum = 8 - 2 * (hashNum - 1);
  fprintf(stdout, "[RAPID] num of hash table changes from %d to %d \n",
          curHashNum, newHashNum);
  curHashNum = newHashNum;
}

int ferret_power(int base, int power){
  int i = 0;
  int result = 1;
  for (;i<power;i++){
    result = result*base;
  }
  return result;
}

void *change_Hash_Num_Cont(void *arg) {
  int hashNum = ferret_power(2,getParaVal(hashContPara));
  fprintf(stdout, "[RAPID] num of hash table changes from %d to %d \n",
          curHashNum, hashNum);
  curHashNum = hashNum;
}

void *change_Probe_Num(void *arg) {
  int probeNum = getParaVal(probePara);
  int newProbeNum;
  newProbeNum = 20 - 2 * (probeNum - 1);
  fprintf(stdout, "[RAPID] num of probing changes from %d to %d \n",
          curProbeNum, newProbeNum);
  curProbeNum = newProbeNum;
}

void *change_Probe_Num_Cont(void *arg) {
  int probeNum = getParaVal(probeContPara);
  fprintf(stdout, "[RAPID] num of hash table changes from %d to %d \n",
          curProbeNum, probeNum);
  curProbeNum = probeNum;
}

// a run function being registered to RSDG
void *change_Iteration_Num(void *arg) {
  int itrNum = getParaVal(iterationPara);
  int newItr;
  newItr = 25 - (itrNum - 1);
  fprintf(stdout, "[RAPID] num of iteration changes from %d to %d\n", curItr,
          newItr);
  curItr = newItr;
}

void *change_Iteration_Num_Cont(void *arg) {
  int iterationNum = getParaVal(iterationContPara);
  fprintf(stdout, "[RAPID] num of iteration changes from %d to %d \n", curItr,
          iterationNum);
  curItr = iterationNum;
}

void setupMission() {
  // init a rsdgmission
  ferretMission = newRAPIDMission(config_file, 1);
  // init a rsdgPara
  iterationPara = newRAPIDPara();
  hashPara = newRAPIDPara();
  probePara = newRAPIDPara();
  iterationContPara = newRAPIDPara();
  hashContPara = newRAPIDPara();
  probeContPara = newRAPIDPara();
  // register discrete service
  if (CONT == 0) {
    for (int i = 0; i < 4; i++) {
      char *node_name = malloc(10 * sizeof(char));
      sprintf(node_name, "hash%d\0", i + 1);
      regService(ferretMission, "hashNum", node_name, &change_Hash_Num, 1,
                 hashPara, i + 1);
    }

    // register the probing service
    for (int i = 0; i < 10; i++) {
      char *node_name = malloc(10 * sizeof(char));
      sprintf(node_name, "probe%d\0", i + 1);
      regService(ferretMission, "probeNum", node_name, &change_Probe_Num, 1,
                 probePara, i + 1);
    }

    // register the iteration service
    for (int i = 0; i < 25; i++) {
      char *node_name = malloc(10 * sizeof(char));
      sprintf(node_name, "itr%d\0", i + 1);
      regService(ferretMission, "iterationNum", node_name,
                 &change_Iteration_Num, 1, iterationPara, i + 1);
    }
  } else { // Continuous Services
    regContService(ferretMission, "hashNum", "hash", &change_Hash_Num_Cont,
                   hashContPara);
    regContService(ferretMission, "probeNum", "probe", &change_Probe_Num_Cont,
                   probeContPara);
    regContService(ferretMission, "iterationNum", "itr",
                   &change_Iteration_Num_Cont, iterationContPara);
  }
  addConstraint(ferretMission, "iterationNum", 1);
  addConstraint(ferretMission, "probeNum", 1);
  addConstraint(ferretMission, "hashNum", 1);
  fprintf(stdout, "[RAPID], Mission Setup Done\n");
}
