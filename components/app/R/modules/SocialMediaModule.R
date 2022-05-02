##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## 
##
##


##=================================================================================
##=============================== SHINY MODULES ===================================
##=================================================================================

SocialMediaModuleUI <- function(id) {
  ns <- shiny::NS(id)
  uiOutput(ns("modal"))
}

SocialMediaModule <- function(id, r.show = reactive(0))
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns
    
    msg <- "Hi. I always thought omics analysis was so difficult, but now I am using BigOmics Playground to analyze my own omics data. No coding required. It's so easy and fun! You should really try it! It's open source and there is even a free version. Go and visit BigOmics at www.bigomics.ch\n\n"
    
    urls <- GetSocialMediaSiteLinks_WithShareLinks(
      url = "http://www.bigomics.ch",
      title = "Analyze omics data yourself! ",
      desc = msg
    )
      
    social.buttons <- fillRow(
      actionButton(ns("twitter"),"", icon=icon("twitter"), class="btn-social"),
      actionButton(ns("linkedin"),"", icon=icon("linkedin"), class="btn-social"),                
      actionButton(ns("facebook"),"", icon=icon("facebook"), class="btn-social"),
      actionButton(ns("whatsapp"),"", icon=icon("whatsapp"), class="btn-social"),                
      actionButton(ns("email"),"", icon=icon("envelope"), class="btn-social"),                
      actionButton(ns("telegram"),"", icon=icon("telegram"), class="btn-social"),                
      actionButton(ns("reddit"),"", icon=icon("reddit"), class="btn-social"),
      actionButton(ns("pinterest"),"", icon=icon("pinterest"), class="btn-social"),                
      actionButton(ns("yahoo"),"", icon=icon("yahoo"), class="btn-social"),
      actionButton(ns("skype"),"", icon=icon("skype"), class="btn-social"),
      actionButton(ns("xing"),"", icon=icon("xing"), class="btn-social")                
      ##tags$a(href=urls["line.me"], NULL, icon("line"), style="font-size:30px;", target="_blank")        
    )
    
    modalButton2 <- function(label, icon=NULL) {
      tags$button(type = "button", class = "btn btn-primary", `data-dismiss` = "modal", 
                  `data-bs-dismiss` = "modal", icon, label)
    }
    
    submit.buttons <- fillRow(
      height = 30,
      flex = c(NA,1,NA),
      actionButton(ns("later"),"Maybe later...", class="btn-warn"),
      br(),
      actionButton(ns("sure"),"Yes, for sure! I just did!", class="btn-primary")
      ##modalButton2(""Yes, sure! I just did!"")
    )  
    
    output$modal <- renderUI({

      do.show <- r.show()
      message("do.show = ", do.show)
      if(do.show==0) return(NULL)
      
      showModal( modalDialog(
        div(HTML("<center><h2>Sorry, time's up mate!</h2></center>"),style="margin-top:0px;"),
        br(),
        HTML("Your FREE session has expired. 
         Did you enjoy using BigOmics Playground? You can extend
         your FREE session by referring BigOmics to your friends!"),
        br(),br(),
        social.buttons,
        br(), br(), br(),
        footer = submit.buttons,
        size = "m",
        easyClose = FALSE
      ))

      shinyjs::disable("sure")

    })

    ## counter for number of referred
    nreferred = 0
    
    ## observe social buttons ------------------------------------

    observeEvent( input$twitter, {
      nreferred <<- nreferred + 1
      browseURL( urls["twitter"] )
      Sys.sleep(5); shinyjs::enable("sure")
    })

    observeEvent( input$linkedin, {
      nreferred <<- nreferred + 1
      browseURL( urls["linkedin"] )
      Sys.sleep(5); shinyjs::enable("sure")      
    })

    observeEvent( input$facebook, {
      nreferred <<- nreferred + 1
      browseURL( urls["facebook"] )
      Sys.sleep(5); shinyjs::enable("sure")            
    })

    observeEvent( input$whatsapp, {
      nreferred <<- nreferred + 1
      browseURL( urls["whatsapp"] )
      Sys.sleep(5); shinyjs::enable("sure")                  
    })

    observeEvent( input$email, {
      nreferred <<- nreferred + 1
      browseURL( urls["email"] )
      Sys.sleep(5); shinyjs::enable("sure")                        
    })

    observeEvent( input$telegram, {
      nreferred <<- nreferred + 1
      browseURL( urls["telegram.me"] )
      Sys.sleep(5); shinyjs::enable("sure")      
    })

    observeEvent( input$reddit, {
      nreferred <<- nreferred + 1
      browseURL( urls["reddit"] )
      Sys.sleep(5); shinyjs::enable("sure")            
    })

    observeEvent( input$pinterest, {
      nreferred <<- nreferred + 1
      browseURL( urls["pinterest"] )
      Sys.sleep(5); shinyjs::enable("sure")      
    })

    observeEvent( input$yahoo, {
      nreferred <<- nreferred + 1
      browseURL( urls["yahoo"] )
      Sys.sleep(5); shinyjs::enable("sure")            
    })

    observeEvent( input$skype, {
      nreferred <<- nreferred + 1
      browseURL( urls["skype"] )
      Sys.sleep(5); shinyjs::enable("sure")            
    })

    observeEvent( input$xing, {
      nreferred <<- nreferred + 1
      browseURL( urls["xing"] )
      Sys.sleep(5); shinyjs::enable("sure")            
    })
    

    ## observe submit buttons ------------------------------------
    rv <- reactiveValues(
        success = 1
    )

    observeEvent( input$later, {      
      rv$success <- 0
      removeModal()
    })
    
    observeEvent( input$sure, {
      if(nreferred==0) return()
      rv$success <- rv$success + 1
      removeModal()
    })

    ## return object --------------------------------
    list(
        success = reactive({ rv$success })
    )
    
  })  ## moduleServer
}


if(FALSE) {
    
    require(shiny)
    shinyApp(
        ui = fluidPage(                   
            actionButton("show","show")
            ##SocialMediaModuleUI("social")
        ),
        server = function(input, output, session) {
            SocialMediaModule(
                "social",
                r.show = reactive(input$show)
            )
        }
    )
    

}

##=================================================================================
## Code apated from https://github.com/bradvin/social-share-urls
##=================================================================================


GetSocialMediaSites_NiceNames <- function() {
  nice_names <- c(
    'add.this'='AddThis',
    'blogger'='Blogger',
    'buffer'='Buffer',
    'diaspora'='Diaspora',
    'douban'='Douban',
    'email'='EMail',
    'evernote'='EverNote',
    'getpocket'='Pocket',
    'facebook'='FaceBook',
    'flattr'='Flattr',
    'flipboard'='FlipBoard',
    'google.bookmarks'='GoogleBookmarks',
    'instapaper'='InstaPaper',
    'line.me'='Line.me',
    'linkedin'='LinkedIn',
    'livejournal'='LiveJournal',
    'gmail'='GMail',
    'hacker.news'='HackerNews',
    'ok.ru'='OK.ru',
    'pinterest'='Pinterest',
    'qzone'='QZone',
    'reddit'='Reddit',
    'renren'='RenRen',
    'skype'='Skype',
    'sms'='SMS',
    'surfingbird.ru'='SurfingBird.ru',
    'telegram.me'='Telegram.me',
    'threema'='Threema',
    'tumblr'='Tumblr',
    'twitter'='Twitter',
    'vk'='VK',
    'weibo'='Weibo',
    'whatsapp'='WhatsApp',
    'xing'='Xing',
    'yahoo'='Yahoo'
  )
  return(nice_names)  
}

GetSocialMediaSites_WithShareLinks_OrderedByPopularity <- function()
{
  sites.ordered <- c(
    'google.bookmarks',
    'facebook',
    'reddit',
    'whatsapp',
    'twitter',
    'linkedin',
    'tumblr',
    'pinterest',
    'blogger',
    'livejournal',
    'evernote',
    'add.this',
    'getpocket',
    'hacker.news',
    'buffer',
    'flipboard',
    'instapaper',
    'surfingbird.ru',
    'flattr',
    'diaspora',
    'qzone',
    'vk',
    'weibo',
    'ok.ru',
    'douban',
    'xing',
    'renren',
    'threema',
    'sms',
    'line.me',
    'skype',
    'telegram.me',
    'email',
    'gmail',
    'yahoo'
  )
  return(sites.ordered)
}

GetSocialMediaSites_WithShareLinks_OrderedByAlphabet <- function() {
  sites <- GetSocialMediaSites_WithShareLinks_OrderedByPopularity()
  return(sort(sites))
}
    
GetSocialMediaSiteLinks_WithShareLinks <- function(
   url="",
   title="",
   image="",
   desc="",
   appid="",
   redirecturl="",
   via="",
   hash_tags="",
   provider="",
   language="",
   user_id="",
   category="",
   phone_number="",
   email_address="",
   cc_email_address="",
   bcc_email_address=""
   )
{      
  
  text <- title        
  text <- paste(text, desc)
  
  urls <- c(
    'add.this' = paste0('http://www.addthis.com/bookmark.php?url=',url),
    'blogger' = paste0('https://www.blogger.com/blog-this.g?u=',url,'&n=',title,'&t=',desc),
    'buffer' = paste0('https://buffer.com/add?text=',text,'&url=',url),
    'diaspora' = paste0('https://share.diasporafoundation.org/?title=','title','&url=',url),
    'douban' = paste0('http://www.douban.com/recommend/?url=',url,'&title=',text),
    'email' = paste0('mailto:',email_address,'?subject=',title,'&body=',desc),
    'evernote' = paste0('https://www.evernote.com/clip.action?url=',url,'&title=',text),
    'getpocket' = paste0('https://getpocket.com/edit?url=',url),
    'facebook' = paste0('http://www.facebook.com/sharer.php?u=',url),
    'flattr' = paste0('https://flattr.com/submit/auto?user_id=',user_id,'&url=',url,'&title=',title,'&description=',text,'&language=',language,'&tags=',hash_tags,'&hidden=HIDDEN&category=',category),
    'flipboard' = paste0('https://share.flipboard.com/bookmarklet/popout?v=2&title=',text,'&url=',url), 
    'gmail' = paste0('https://mail.google.com/mail/?view=cm&to=',email_address,'&su=',title,'&body=',url,'&bcc=',bcc_email_address,'&cc=',cc_email_address),
    'google.bookmarks' = paste0('https://www.google.com/bookmarks/mark?op=edit&bkmk=',url,'&title=',title,'&annotation=',text,'&labels=',hash_tags,''),
    'instapaper' = paste0('http://www.instapaper.com/edit?url=',url,'&title=',title,'&description=',desc),
    'line.me' = paste0('https://lineit.line.me/share/ui?url=',url,'&text=',text),
    'linkedin' = paste0('https://www.linkedin.com/sharing/share-offsite/?url=',url),
    'livejournal' = paste0('http://www.livejournal.com/update.bml?subject=',text,'&event=',url),
    'hacker.news' = paste0('https://news.ycombinator.com/submitlink?u=',url,'&t=',title),
    'ok.ru' = paste0('https://connect.ok.ru/dk?st.cmd=WidgetSharePreview&st.shareUrl=',url),
    'pinterest' = paste0('http://pinterest.com/pin/create/button/?url=',url),
    'qzone' = paste0('http://sns.qzone.qq.com/cgi-bin/qzshare/cgi_qzshare_onekey?url=',url),
    'reddit' = paste0('https://reddit.com/submit?url=',url,'&title=',title),
    'renren' = paste0('http://widget.renren.com/dialog/share?resourceUrl=',url,'&srcUrl=',url,'&title=',
                      text,'&description=',desc),
    'skype' = paste0('https://web.skype.com/share?url=',url,'&text=',text),
    'sms' = paste0('sms:',phone_number,'?body=',text),
    'surfingbird.ru' = paste0('http://surfingbird.ru/share?url=',url,'&description=',desc,'&screenshot=',
                              image,'&title=',title),
    'telegram.me' = paste0('https://t.me/share/url?url=',url,'&text=',text,'&to=',phone_number),
    'threema' = paste0('threema://compose?text=',text,'&id=',user_id),
    'tumblr' = paste0('https://www.tumblr.com/widgets/share/tool?canonicalUrl=',url,'&title=',title,
                      '&caption=',desc,'&tags=',hash_tags),
    'twitter' = paste0('https://twitter.com/intent/tweet?url=',url,'&text=',text,'&via=',
                       via,'&hashtags=',hash_tags),
    'vk' = paste0('http://vk.com/share.php?url=',url,'&title=',title,'&comment=',desc),
    'weibo' = paste0('http://service.weibo.com/share/share.php?url=',url,'&appkey=&title=',title,'&pic=&ralateUid='),
    'whatsapp' = paste0('https://api.whatsapp.com/send?text=',text,'%20',url),
    'xing' = paste0('https://www.xing.com/spi/shares/new?url=',url),
    'yahoo' = paste0('http://compose.mail.yahoo.com/?to=','email_address','&subject=',title,'&body=',text)
  )
  return(urls)
}


if(FALSE) {
  
  socialmediasites = GetSocialMediaSites_WithShareLinks_OrderedByAlphabet()
  head(socialmediasites)

  socialmediasites = GetSocialMediaSites_WithShareLinks_OrderedByPopularity()
  head(socialmediasites)

  socialmediaurls <- GetSocialMediaSiteLinks_WithShareLinks(
    url = 'http://www.earthfluent.com/',
    title = 'EarthFluent'
  )
  head(socialmediaurls)
  
  for (socialmediasite in socialmediasites) {
    cat(socialmediasite,":", socialmediaurls[socialmediasite],"\n")
  }

}
