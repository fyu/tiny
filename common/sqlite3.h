#ifndef FURRY_COMMON_SQLIT3_H
#define FURRY_COMMON_SQLIT3_H

#include "furry/3rdparty/sqlite3/sqlite3.h"

#define SQLITE3_BEGIN_TRANSACTION_OR_DIE(db)                            \
  do {                                                                  \
    char *errmsg = 0;                                                   \
    int rc;                                                             \
    rc = sqlite3_exec((db), "BEGIN TRANSACTION", NULL, NULL, &errmsg);  \
    if (rc != SQLITE_OK) {                                              \
      LOG(FATAL) << "Failed to begin transation. "                      \
                 << "SQL error: " << errmsg;                            \
      sqlite3_free(errmsg);                                             \
    }                                                                   \
  } while(0)                                                            \

#define SQLITE3_END_TRANSACTION_OR_DIE(db)                            \
  do {                                                                \
    char *errmsg = 0;                                                 \
    int rc;                                                           \
    rc = sqlite3_exec((db), "END TRANSACTION", NULL, NULL, &errmsg);  \
    if (rc != SQLITE_OK) {                                            \
      LOG(FATAL) << "Failed to end transation. "                      \
                 << "SQL error: " << errmsg;                          \
      sqlite3_free(errmsg);                                           \
    }                                                                 \
  } while(0)                                                          \

#define SQLITE3_EXEC_OR_DIR(db, sql)                          \
  do {                                                        \
    char *errmsg = 0;                                         \
    int rc = sqlite3_exec((db), (sql), nullptr, 0, &errmsg);  \
    if (rc != SQLITE_OK) {                                    \
      LOG(FATAL) << "Failed to exec transation. "             \
                 << "SQL error: " << errmsg;                  \
      sqlite3_free(errmsg);                                   \
    }                                                         \
  } while(0)


#endif // FURRY_COMMON_SQLIT3_H
